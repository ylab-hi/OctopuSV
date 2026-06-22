"""SVCF sample/caller subset engine.

subset = evidence-column operation

It is intentionally different from:
    filter = row-level SV attribute filtering
    query  = row-level target-based retrieval

OctopuSV SVCF has two important layouts:

1. sample mode
   Header:
       #CHROM ... FORMAT sample1 sample2 sample3

   Meaning:
       FORMAT trailing columns are sample/input columns.

2. caller mode
   Header:
       #CHROM ... FORMAT SAMPLE

   But each record may contain multiple trailing evidence blocks:
       FORMAT SAMPLE
       [caller1 evidence] [caller2 evidence] ...

   Meaning:
       INFO/SOURCES defines the caller order and maps to the trailing evidence
       blocks. This is SVCF-specific and must not be treated as ordinary VCF.

Core rules:
    - sample mode uses #CHROM sample column names as truth.
    - caller mode uses INFO/SOURCES order as truth.
    - SOURCES and SOURCE_IDS are updated after subset.
    - SUPPORT is not recomputed.
    - By default, records unsupported by the selected sample/caller set are
      dropped.
    - --keep-empty keeps empty shell records.
"""

from __future__ import annotations

import json
from collections import Counter
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

from octopusv.filtering.svcf_filter import parse_info


LEGAL_MODES = {"auto", "sample", "caller"}
MISSING_VALUES = {None, "", "."}


@dataclass
class SubsetConfig:
    input_file: str
    output_file: Optional[str] = None
    dry_run: bool = False

    mode: str = "auto"
    selected_samples: list[str] = field(default_factory=list)
    selected_callers: list[str] = field(default_factory=list)

    keep_empty: bool = False
    audit_header: bool = True

    def __post_init__(self) -> None:
        self.mode = str(self.mode).lower()

        if self.mode not in LEGAL_MODES:
            raise ValueError(
                f"Invalid --mode '{self.mode}'. "
                f"Allowed values: {sorted(LEGAL_MODES)}."
            )

        self.selected_samples = unique_preserve(self.selected_samples)
        self.selected_callers = unique_preserve(self.selected_callers)

        if self.selected_samples and self.selected_callers:
            raise ValueError("Use either --sample or --caller, not both.")

        if not self.selected_samples and not self.selected_callers:
            raise ValueError("At least one --sample or --caller is required.")

        if not self.dry_run and self.output_file is None:
            raise ValueError("output_file is required unless dry_run=True.")


@dataclass
class HeaderInfo:
    meta_lines: list[str]
    chrom_line: str
    chrom_fields: list[str]
    has_multi_mode_tag: bool
    detected_mode: str
    mode_warnings: list[str]
    available_samples: list[str]


def unique_preserve(values: list[str]) -> list[str]:
    seen = set()
    result = []

    for value in values:
        value = str(value).strip()
        if not value:
            continue
        if value in seen:
            continue
        seen.add(value)
        result.append(value)

    return result


def expand_csv_values(values: Optional[list[str]]) -> list[str]:
    """Support both repeated options and comma-separated values."""
    if not values:
        return []

    expanded: list[str] = []
    for value in values:
        for item in str(value).split(","):
            item = item.strip()
            if item:
                expanded.append(item)

    return unique_preserve(expanded)


def read_name_file(path: Optional[Path | str]) -> list[str]:
    if path is None:
        return []

    names: list[str] = []
    with Path(path).open() as handle:
        for line in handle:
            value = line.strip()
            if not value or value.startswith("#"):
                continue
            names.append(value)

    return unique_preserve(names)


def split_info_list(value) -> list[str]:
    if value in MISSING_VALUES or value is True:
        return []

    items = []
    for item in str(value).split(","):
        item = item.strip()
        if item and item != ".":
            items.append(item)

    return items


def update_info_string(info_string: str, updates: dict[str, str]) -> str:
    """Update INFO while preserving existing field order as much as possible."""
    if info_string in MISSING_VALUES:
        items: list[str] = []
    else:
        items = [item for item in str(info_string).split(";") if item]

    key_to_index: dict[str, int] = {}
    for index, item in enumerate(items):
        key = item.split("=", 1)[0]
        key_to_index[key] = index

    for key, value in updates.items():
        value = "." if value in MISSING_VALUES else str(value)
        new_item = f"{key}={value}"

        if key in key_to_index:
            items[key_to_index[key]] = new_item
        else:
            items.append(new_item)

    return ";".join(items) if items else "."


def evidence_is_supported(evidence: str) -> bool:
    """Return whether a sample/caller evidence field looks supportive.

    SVCF FORMAT starts with GT in normal OctopuSV output.
    Unsupported placeholders are usually 0/0 or missing-like values.
    """
    if evidence in MISSING_VALUES:
        return False

    text = str(evidence).strip()
    if not text or text == ".":
        return False

    gt = text.split(":", 1)[0]
    if gt in {"0/0", "0|0", "./.", ".|.", "."}:
        return False

    return True


def make_empty_evidence(format_field: str) -> str:
    """Create a minimal 0/0 placeholder matching FORMAT length."""
    if format_field in MISSING_VALUES:
        return "0/0"

    keys = str(format_field).split(":")
    if not keys:
        return "0/0"

    return ":".join(["0/0"] + ["."] * (len(keys) - 1))


class SVCFSubset:
    """Subset SVCF sample/caller evidence columns."""

    def __init__(self, config: SubsetConfig):
        self.config = config

        self.input_records = 0
        self.output_records = 0
        self.dropped_empty_records = 0
        self.empty_records_retained = 0

        self.detected_mode: Optional[str] = None
        self.selected_type: Optional[str] = None

        self.available_samples: list[str] = []
        self.retained_samples: list[str] = []

        self.available_callers_counter = Counter()
        self.retained_callers_counter = Counter()

        self.warning_counts = Counter()
        self.mode_warnings: list[str] = []

        self.header_info: Optional[HeaderInfo] = None

    def run(self) -> dict:
        input_path = Path(self.config.input_file)

        if not input_path.exists():
            raise FileNotFoundError(f"Input file not found: {input_path}")

        output_handle = None
        if not self.config.dry_run:
            output_handle = Path(self.config.output_file).open("w")

        try:
            with input_path.open() as handle:
                for line in handle:
                    if line.startswith("##"):
                        self._collect_meta_header_line(line)
                        continue

                    if line.startswith("#CHROM"):
                        self._finalize_header(line)

                        if output_handle is not None:
                            for header_line in self._render_output_header_lines():
                                output_handle.write(header_line)

                        continue

                    if not line.strip():
                        continue

                    if self.header_info is None:
                        raise ValueError("Malformed SVCF: data record found before #CHROM header.")

                    self.input_records += 1

                    output_line = self._process_record(line)
                    if output_line is None:
                        continue

                    self.output_records += 1

                    if output_handle is not None:
                        output_handle.write(output_line)

        finally:
            if output_handle is not None:
                output_handle.close()

        if self.header_info is None:
            raise ValueError("Malformed SVCF: missing #CHROM header line.")

        return self.summary()

    def _collect_meta_header_line(self, line: str) -> None:
        if self.header_info is None:
            if not hasattr(self, "_meta_lines"):
                self._meta_lines = []

            # Avoid accumulating stale subset audit lines after repeated subset.
            if line.startswith("##OctopuSV_subset_"):
                return

            self._meta_lines.append(line)

    def _finalize_header(self, chrom_line: str) -> None:
        meta_lines = getattr(self, "_meta_lines", [])
        chrom_fields = chrom_line.rstrip("\n").split("\t")

        if len(chrom_fields) < 9:
            raise ValueError("Malformed SVCF header: #CHROM line has fewer than 9 columns.")

        has_multi_mode_tag = any(
            line.strip() == "##OctopuSV_mode=multi" for line in meta_lines
        )

        available_samples = chrom_fields[9:]
        detected_mode, mode_warnings = self._detect_mode(
            requested_mode=self.config.mode,
            has_multi_mode_tag=has_multi_mode_tag,
            available_samples=available_samples,
        )

        self.detected_mode = detected_mode
        self.mode_warnings = mode_warnings
        for warning in mode_warnings:
            self.warning_counts[warning] += 1

        self.header_info = HeaderInfo(
            meta_lines=meta_lines,
            chrom_line=chrom_line,
            chrom_fields=chrom_fields,
            has_multi_mode_tag=has_multi_mode_tag,
            detected_mode=detected_mode,
            mode_warnings=mode_warnings,
            available_samples=available_samples,
        )

        self._validate_selection_against_mode()

    def _detect_mode(
        self,
        requested_mode: str,
        has_multi_mode_tag: bool,
        available_samples: list[str],
    ) -> tuple[str, list[str]]:
        """Detect sample/caller mode using positive header signatures."""
        warnings: list[str] = []

        trailing_count = len(available_samples)

        if trailing_count < 1:
            raise ValueError("Malformed SVCF header: no FORMAT trailing column found.")

        # Explicit mode request.
        if requested_mode in {"sample", "caller"}:
            if requested_mode == "caller" and trailing_count != 1:
                raise ValueError(
                    "--mode caller requires a caller-mode header with exactly one "
                    "FORMAT trailing column."
                )

            if requested_mode == "caller" and has_multi_mode_tag:
                warnings.append("forced_caller_mode_despite_multi_mode_tag")

            if requested_mode == "sample" and trailing_count == 1 and available_samples[0] == "SAMPLE":
                warnings.append("forced_sample_mode_with_single_SAMPLE_header")

            return requested_mode, warnings

        # Auto mode.
        if has_multi_mode_tag:
            if trailing_count == 1 and available_samples[0] == "SAMPLE":
                raise ValueError(
                    "Header conflict: ##OctopuSV_mode=multi is present but #CHROM "
                    "has only a single SAMPLE column. Use --mode explicitly after "
                    "checking the file."
                )
            return "sample", warnings

        if trailing_count == 1 and available_samples[0] == "SAMPLE":
            return "caller", warnings

        # Multi-column header is a strong sample-mode signature even if the
        # OctopuSV mode tag was lost by a previous operation.
        if trailing_count > 1:
            warnings.append("mode_tag_missing_but_multicolumn_header_inferred_sample_mode")
            return "sample", warnings

        # Single named column other than SAMPLE behaves like sample-mode.
        warnings.append("mode_tag_missing_single_named_sample_inferred_sample_mode")
        return "sample", warnings

    def _validate_selection_against_mode(self) -> None:
        assert self.header_info is not None

        if self.detected_mode == "sample":
            if self.config.selected_callers:
                raise ValueError(
                    "This file is detected as sample-mode SVCF. Use --sample, not --caller."
                )

            self.selected_type = "sample"
            available = self.header_info.available_samples
            available_set = set(available)

            missing = [s for s in self.config.selected_samples if s not in available_set]
            if missing:
                raise ValueError(
                    "Requested sample(s) not found in #CHROM header: "
                    + ",".join(missing)
                )

            self.available_samples = available
            self.retained_samples = [
                sample for sample in available if sample in set(self.config.selected_samples)
            ]

            if not self.retained_samples:
                raise ValueError("No selected samples remain after header matching.")

        elif self.detected_mode == "caller":
            if self.config.selected_samples:
                raise ValueError(
                    "This file is detected as caller-mode SVCF. Use --caller, not --sample."
                )

            self.selected_type = "caller"
            if not self.config.selected_callers:
                raise ValueError("--caller is required for caller-mode SVCF.")

        else:
            raise ValueError(f"Internal error: unknown detected mode {self.detected_mode}")

    def _render_output_header_lines(self) -> list[str]:
        assert self.header_info is not None

        lines: list[str] = []
        lines.extend(self.header_info.meta_lines)

        if self.config.audit_header:
            selected = (
                self.retained_samples
                if self.detected_mode == "sample"
                else self.config.selected_callers
            )
            lines.append(f"##OctopuSV_subset_mode={self.detected_mode}\n")
            lines.append(f"##OctopuSV_subset_selected_type={self.selected_type}\n")
            lines.append(f"##OctopuSV_subset_selected={','.join(selected)}\n")
            lines.append(f"##OctopuSV_subset_keep_empty={str(self.config.keep_empty).lower()}\n")

        fields = self.header_info.chrom_fields

        if self.detected_mode == "sample":
            output_fields = fields[:9] + self.retained_samples
        else:
            # Caller mode keeps the single header sample column unchanged.
            output_fields = fields

        return lines + ["\t".join(output_fields) + "\n"]

    def _process_record(self, line: str) -> Optional[str]:
        if self.detected_mode == "sample":
            return self._process_sample_mode_record(line)

        if self.detected_mode == "caller":
            return self._process_caller_mode_record(line)

        raise ValueError("Internal error: mode was not detected before record processing.")

    def _process_sample_mode_record(self, line: str) -> Optional[str]:
        assert self.header_info is not None

        fields = line.rstrip("\n").split("\t")
        if len(fields) < 9:
            self.warning_counts["malformed_record"] += 1
            return None

        fixed = fields[:9]
        info_string = fields[7]
        format_field = fields[8]
        evidence_fields = fields[9:]

        available_samples = self.header_info.available_samples
        selected_set = set(self.retained_samples)

        if len(evidence_fields) != len(available_samples):
            self.warning_counts["record_sample_column_count_mismatch"] += 1

        padded_evidence = evidence_fields[:]
        while len(padded_evidence) < len(available_samples):
            padded_evidence.append(make_empty_evidence(format_field))

        padded_evidence = padded_evidence[: len(available_samples)]

        selected_indices = [
            index
            for index, sample in enumerate(available_samples)
            if sample in selected_set
        ]

        selected_evidence = [padded_evidence[index] for index in selected_indices]

        retained_sources = []
        for sample, evidence in zip(self.retained_samples, selected_evidence):
            if evidence_is_supported(evidence):
                retained_sources.append(sample)

        if not retained_sources and not self.config.keep_empty:
            self.dropped_empty_records += 1
            return None

        if not retained_sources and self.config.keep_empty:
            self.empty_records_retained += 1

        info = parse_info(info_string)
        source_ids_update = self._subset_source_ids_for_sample_mode(
            info=info,
            retained_sources=retained_sources,
        )

        updates = {
            "SOURCES": ",".join(retained_sources) if retained_sources else ".",
        }

        if "SOURCE_IDS" in info:
            updates["SOURCE_IDS"] = source_ids_update

        fixed[7] = update_info_string(info_string, updates)

        return "\t".join(fixed + selected_evidence) + "\n"

    def _subset_source_ids_for_sample_mode(
        self,
        info: dict,
        retained_sources: list[str],
    ) -> str:
        sources = split_info_list(info.get("SOURCES"))
        source_ids = split_info_list(info.get("SOURCE_IDS"))

        if "SOURCE_IDS" not in info:
            return "."

        if not retained_sources:
            return "."

        if len(sources) != len(source_ids):
            self.warning_counts["source_ids_count_mismatch"] += 1
            return "."

        source_to_id = dict(zip(sources, source_ids))

        retained_ids = []
        missing_ids = False
        for source in retained_sources:
            if source not in source_to_id:
                missing_ids = True
                continue
            retained_ids.append(source_to_id[source])

        if missing_ids:
            self.warning_counts["selected_source_missing_source_id"] += 1

        return ",".join(retained_ids) if retained_ids else "."

    def _process_caller_mode_record(self, line: str) -> Optional[str]:
        fields = line.rstrip("\n").split("\t")
        if len(fields) < 9:
            self.warning_counts["malformed_record"] += 1
            return None

        fixed = fields[:9]
        info_string = fields[7]
        format_field = fields[8]
        evidence_fields = fields[9:]

        info = parse_info(info_string)
        sources = split_info_list(info.get("SOURCES"))

        for source in sources:
            self.available_callers_counter[source] += 1

        selected_set = set(self.config.selected_callers)

        if not sources:
            self.warning_counts["sources_missing_in_caller_mode"] += 1
            if not self.config.keep_empty:
                self.dropped_empty_records += 1
                return None

            selected_evidence = [
                make_empty_evidence(format_field) for _ in self.config.selected_callers
            ]
            self.empty_records_retained += 1

            updates = {"SOURCES": "."}
            if "SOURCE_IDS" in info:
                updates["SOURCE_IDS"] = "."

            fixed[7] = update_info_string(info_string, updates)
            return "\t".join(fixed + selected_evidence) + "\n"

        if len(evidence_fields) != len(sources):
            self.warning_counts["source_evidence_count_mismatch"] += 1
            if not self.config.keep_empty:
                self.dropped_empty_records += 1
                return None

            selected_evidence = [
                make_empty_evidence(format_field) for _ in self.config.selected_callers
            ]
            self.empty_records_retained += 1

            updates = {"SOURCES": "."}
            if "SOURCE_IDS" in info:
                updates["SOURCE_IDS"] = "."

            fixed[7] = update_info_string(info_string, updates)
            return "\t".join(fixed + selected_evidence) + "\n"

        retained_indices = [
            index
            for index, source in enumerate(sources)
            if source in selected_set
        ]

        retained_sources = [sources[index] for index in retained_indices]
        selected_evidence = [evidence_fields[index] for index in retained_indices]

        for source in retained_sources:
            self.retained_callers_counter[source] += 1

        if not retained_sources and not self.config.keep_empty:
            self.dropped_empty_records += 1
            return None

        if not retained_sources and self.config.keep_empty:
            self.empty_records_retained += 1
            selected_evidence = [
                make_empty_evidence(format_field) for _ in self.config.selected_callers
            ]

        for source, evidence in zip(retained_sources, selected_evidence):
            if not evidence_is_supported(evidence):
                self.warning_counts["selected_source_has_empty_evidence"] += 1

        source_ids_update = self._subset_source_ids_for_caller_mode(
            info=info,
            sources=sources,
            retained_indices=retained_indices,
        )

        updates = {
            "SOURCES": ",".join(retained_sources) if retained_sources else ".",
        }

        if "SOURCE_IDS" in info:
            updates["SOURCE_IDS"] = source_ids_update

        fixed[7] = update_info_string(info_string, updates)

        return "\t".join(fixed + selected_evidence) + "\n"

    def _subset_source_ids_for_caller_mode(
        self,
        info: dict,
        sources: list[str],
        retained_indices: list[int],
    ) -> str:
        source_ids = split_info_list(info.get("SOURCE_IDS"))

        if "SOURCE_IDS" not in info:
            return "."

        if not retained_indices:
            return "."

        if len(sources) != len(source_ids):
            self.warning_counts["source_ids_count_mismatch"] += 1
            return "."

        retained_ids = [source_ids[index] for index in retained_indices]
        return ",".join(retained_ids) if retained_ids else "."

    def summary(self) -> dict:
        selected = (
            self.retained_samples
            if self.detected_mode == "sample"
            else self.config.selected_callers
        )

        return {
            "tool": "octopusv subset",
            "mode": "evidence_column_subset",
            "input": self.config.input_file,
            "output": None if self.config.dry_run else self.config.output_file,
            "dry_run": self.config.dry_run,
            "requested_mode": self.config.mode,
            "detected_mode": self.detected_mode,
            "selected_type": self.selected_type,
            "selected": selected,
            "keep_empty": self.config.keep_empty,
            "input_records": self.input_records,
            "output_records": self.output_records,
            "dropped_empty_records": self.dropped_empty_records,
            "empty_records_retained": self.empty_records_retained,
            "preserves_header": True,
            "updates_sources": True,
            "updates_source_ids": True,
            "recomputes_support": False,
            "available_samples": self.available_samples,
            "retained_samples": self.retained_samples,
            "observed_callers": sorted(self.available_callers_counter),
            "retained_callers": sorted(self.retained_callers_counter),
            "warning_counts": dict(self.warning_counts),
            "mode_warnings": self.mode_warnings,
        }

    @staticmethod
    def summary_to_json(summary: dict) -> str:
        return json.dumps(summary, indent=2)

    @staticmethod
    def summary_to_text(summary: dict) -> str:
        lines = [
            "SVCF subset summary",
            f"Input: {summary['input']}",
            f"Output: {summary['output']}",
            f"Dry run: {summary['dry_run']}",
            f"Requested mode: {summary['requested_mode']}",
            f"Detected mode: {summary['detected_mode']}",
            f"Selected type: {summary['selected_type']}",
            f"Selected: {','.join(summary['selected'])}",
            f"Keep empty: {summary['keep_empty']}",
            f"Input records: {summary['input_records']}",
            f"Output records: {summary['output_records']}",
            f"Dropped empty records: {summary['dropped_empty_records']}",
            f"Empty records retained: {summary['empty_records_retained']}",
            f"Updates SOURCES: {summary['updates_sources']}",
            f"Updates SOURCE_IDS: {summary['updates_source_ids']}",
            f"Recomputes SUPPORT: {summary['recomputes_support']}",
        ]

        if summary["available_samples"]:
            lines.append("")
            lines.append("Available samples:")
            for sample in summary["available_samples"]:
                lines.append(f"  {sample}")

        if summary["retained_samples"]:
            lines.append("")
            lines.append("Retained samples:")
            for sample in summary["retained_samples"]:
                lines.append(f"  {sample}")

        if summary["observed_callers"]:
            lines.append("")
            lines.append("Observed callers:")
            for caller in summary["observed_callers"]:
                lines.append(f"  {caller}")

        if summary["retained_callers"]:
            lines.append("")
            lines.append("Retained callers with output records:")
            for caller in summary["retained_callers"]:
                lines.append(f"  {caller}")

        if summary["mode_warnings"]:
            lines.append("")
            lines.append("Mode warnings:")
            for warning in summary["mode_warnings"]:
                lines.append(f"  {warning}")

        if summary["warning_counts"]:
            lines.append("")
            lines.append("Warnings:")
            for warning, count in sorted(summary["warning_counts"].items()):
                lines.append(f"  {warning}: {count}")

        return "\n".join(lines)