"""SVCF-aware row filtering engine.

This module implements record-level filtering for OctopuSV SVCF files.

Important boundary:
    filter = row-level operation
    subset = column-level operation

Therefore:
    filter --source pbsv
        keeps records supported by pbsv, but preserves all columns.

    subset --caller pbsv
        later should keep only the pbsv caller/evidence column.

This engine preserves all SVCF header lines and writes filtered records
without rewriting INFO/FORMAT/sample/caller columns.
"""

from __future__ import annotations

import json
from collections import Counter
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional


NON_LINEAR_SVTYPES = {"TRA", "BND"}
LEGAL_SOURCE_MODES = {"any", "all"}


def parse_csv_arg(value: Optional[str]) -> set[str]:
    """Parse comma-separated CLI input into a set of non-empty strings."""
    if value is None:
        return set()

    values = set()
    for item in str(value).split(","):
        item = item.strip()
        if item:
            values.add(item)
    return values


def read_list_file(path: Optional[Path | str]) -> set[str]:
    """Read one value per line, ignoring blank lines and comments."""
    if path is None:
        return set()

    values = set()
    with Path(path).open() as handle:
        for line in handle:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            values.add(line)
    return values


def parse_info(info_str: str) -> dict:
    """Parse VCF/SVCF INFO into a dict.

    INFO flags without '=' are stored as True.
    """
    info = {}
    if info_str in ("", "."):
        return info

    for item in info_str.split(";"):
        if not item:
            continue
        if "=" in item:
            key, value = item.split("=", 1)
            info[key] = value
        else:
            info[item] = True
    return info


def split_filter_tokens(filter_value: str) -> set[str]:
    """Split FILTER column by ';'.

    Exact FILTER matching and token matching are different:
        FILTER == "Decoy;NearReferenceGap"
        tokens == {"Decoy", "NearReferenceGap"}
    """
    if filter_value in ("", "."):
        return {filter_value}
    return {token for token in filter_value.split(";") if token}


def safe_float(value) -> Optional[float]:
    """Return float(value), or None for missing/non-numeric values."""
    if value in (None, "", ".", True):
        return None

    try:
        # Some INFO values may be comma-separated. Use the first value.
        first = str(value).split(",", 1)[0]
        return float(first)
    except (TypeError, ValueError):
        return None


def safe_int(value) -> Optional[int]:
    """Return int(value), or None for missing/non-integer values."""
    if value in (None, "", ".", True):
        return None

    try:
        first = str(value).split(",", 1)[0]
        return int(first)
    except (TypeError, ValueError):
        return None


def normalize_standard_contig(contig: str) -> Optional[str]:
    """Normalize standard human contigs.

    Returns:
        "1"-"22", "X", "Y", or "MT" for standard chromosomes.
        None for nonstandard contigs.

    Supports:
        1 / chr1
        X / chrX
        Y / chrY
        M / MT / chrM / chrMT
    """
    if contig is None:
        return None

    c = str(contig).strip()
    if c.lower().startswith("chr"):
        c = c[3:]

    if c in {"M", "MT"}:
        return "MT"

    if c in {"X", "Y"}:
        return c

    if c.isdigit():
        n = int(c)
        if 1 <= n <= 22:
            return str(n)

    return None


def is_standard_contig(contig: str) -> bool:
    return normalize_standard_contig(contig) is not None


def contig_matches(contig: str, targets: set[str]) -> bool:
    """Check exact or standard-chromosome alias match.

    Examples:
        contig_matches("chr1", {"1"}) is True
        contig_matches("1", {"chr1"}) is True
        contig_matches("GL000220.1", {"GL000220.1"}) is True
    """
    if not targets:
        return False

    contig_str = str(contig)
    target_lower = {t.lower() for t in targets}

    if contig_str.lower() in target_lower:
        return True

    contig_norm = normalize_standard_contig(contig_str)
    if contig_norm is None:
        return False

    for target in targets:
        target_norm = normalize_standard_contig(target)
        if target_norm is not None and target_norm == contig_norm:
            return True

    return False


def normalize_svtype_set(values: set[str]) -> set[str]:
    return {v.upper() for v in values}


def normalize_source_set(values: set[str]) -> set[str]:
    return {v.lower() for v in values}


def svcf_event_size(svtype: str, pos: str, info: dict) -> Optional[int]:
    """Return OctopuSV SVCF-defined SV size.

    Rules:
        DEL/DUP/INV/INS:
            use numeric SVLEN if available;
            otherwise use abs(END - POS).

        TRA/BND:
            no linear size; return None.

    In OctopuSV SVCF, INS can have meaningful POS/END span.
    """
    svtype = str(svtype).upper()

    if svtype in NON_LINEAR_SVTYPES:
        return None

    svlen = safe_int(info.get("SVLEN"))
    if svlen is not None:
        return abs(svlen)

    pos_i = safe_int(pos)
    end_i = safe_int(info.get("END"))
    if pos_i is None or end_i is None:
        return None

    return abs(end_i - pos_i)


@dataclass
class FilterConfig:
    """Configuration for SVCF row-level filtering."""

    input_file: str
    output_file: Optional[str] = None
    dry_run: bool = False

    svtypes: set[str] = field(default_factory=set)
    exclude_svtypes: set[str] = field(default_factory=set)

    pass_only: bool = False
    filters: set[str] = field(default_factory=set)
    exclude_filters: set[str] = field(default_factory=set)
    filter_tokens: set[str] = field(default_factory=set)
    exclude_filter_tokens: set[str] = field(default_factory=set)

    min_size: Optional[int] = None
    max_size: Optional[int] = None

    min_support: Optional[int] = None
    max_support: Optional[int] = None

    sources: set[str] = field(default_factory=set)
    exclude_sources: set[str] = field(default_factory=set)
    source_mode: str = "any"
    min_source_count: Optional[int] = None
    max_source_count: Optional[int] = None

    ids: set[str] = field(default_factory=set)
    exclude_ids: set[str] = field(default_factory=set)

    contigs: set[str] = field(default_factory=set)
    exclude_contigs: set[str] = field(default_factory=set)
    standard_contigs: bool = False

    min_qual: Optional[float] = None
    max_qual: Optional[float] = None

    min_af: Optional[float] = None
    max_af: Optional[float] = None

    def __post_init__(self) -> None:
        self.svtypes = normalize_svtype_set(self.svtypes)
        self.exclude_svtypes = normalize_svtype_set(self.exclude_svtypes)

        self.sources = normalize_source_set(self.sources)
        self.exclude_sources = normalize_source_set(self.exclude_sources)

        if self.source_mode not in LEGAL_SOURCE_MODES:
            raise ValueError(
                f"Invalid source mode '{self.source_mode}'. "
                f"Allowed values: {sorted(LEGAL_SOURCE_MODES)}."
            )

        if self.min_size is not None and self.min_size < 0:
            raise ValueError("--min-size must be non-negative.")
        if self.max_size is not None and self.max_size < 0:
            raise ValueError("--max-size must be non-negative.")
        if (
            self.min_size is not None
            and self.max_size is not None
            and self.min_size > self.max_size
        ):
            raise ValueError("--min-size cannot be greater than --max-size.")

        if self.min_support is not None and self.min_support < 0:
            raise ValueError("--min-support must be non-negative.")
        if self.max_support is not None and self.max_support < 0:
            raise ValueError("--max-support must be non-negative.")
        if (
            self.min_support is not None
            and self.max_support is not None
            and self.min_support > self.max_support
        ):
            raise ValueError("--min-support cannot be greater than --max-support.")

    def active_filters(self) -> dict:
        """Return JSON-friendly representation of active filters."""
        return {
            "svtype": sorted(self.svtypes),
            "exclude_svtype": sorted(self.exclude_svtypes),
            "pass_only": self.pass_only,
            "filter": sorted(self.filters),
            "exclude_filter": sorted(self.exclude_filters),
            "filter_token": sorted(self.filter_tokens),
            "exclude_filter_token": sorted(self.exclude_filter_tokens),
            "min_size": self.min_size,
            "max_size": self.max_size,
            "min_support": self.min_support,
            "max_support": self.max_support,
            "source": sorted(self.sources),
            "exclude_source": sorted(self.exclude_sources),
            "source_mode": self.source_mode,
            "min_source_count": self.min_source_count,
            "max_source_count": self.max_source_count,
            "id": sorted(self.ids),
            "exclude_id": sorted(self.exclude_ids),
            "contig": sorted(self.contigs),
            "exclude_contig": sorted(self.exclude_contigs),
            "standard_contigs": self.standard_contigs,
            "min_qual": self.min_qual,
            "max_qual": self.max_qual,
            "min_af": self.min_af,
            "max_af": self.max_af,
        }


class SVCFFilter:
    """Apply SVCF-aware row filters while preserving SVCF structure."""

    def __init__(self, config: FilterConfig):
        self.config = config
        self.header_trailing_columns: list[str] = []

        self.input_records = 0
        self.output_records = 0

        self.excluded_by_reason = Counter()
        self.source_detection_used = Counter()

    def run(self) -> dict:
        """Run filtering and return a summary dict."""
        input_path = Path(self.config.input_file)

        if not input_path.exists():
            raise FileNotFoundError(f"Input file not found: {input_path}")

        output_handle = None
        if not self.config.dry_run:
            if self.config.output_file is None:
                raise ValueError("output_file is required unless dry_run=True.")
            output_handle = Path(self.config.output_file).open("w")

        try:
            with input_path.open() as handle:
                for line in handle:
                    if line.startswith("#"):
                        self._handle_header(line, output_handle)
                        continue

                    if not line.strip():
                        continue

                    self.input_records += 1

                    keep, reason = self._passes_line(line)
                    if keep:
                        self.output_records += 1
                        if output_handle is not None:
                            output_handle.write(line)
                    else:
                        self.excluded_by_reason[reason] += 1
        finally:
            if output_handle is not None:
                output_handle.close()

        return self.summary()

    def _handle_header(self, line: str, output_handle) -> None:
        if line.startswith("#CHROM"):
            columns = line.rstrip("\n").split("\t")
            self.header_trailing_columns = columns[9:] if len(columns) > 9 else []

        if output_handle is not None:
            output_handle.write(line)

    def _passes_line(self, line: str) -> tuple[bool, str]:
        fields = line.rstrip("\n").split("\t")
        if len(fields) < 8:
            return False, "malformed_record"

        chrom = fields[0]
        pos = fields[1]
        sv_id = fields[2]
        qual = fields[5] if len(fields) > 5 else "."
        filter_value = fields[6] if len(fields) > 6 else "."
        info = parse_info(fields[7])

        svtype = str(info.get("SVTYPE", "Unknown")).upper()

        # ID filters
        if self.config.ids and sv_id not in self.config.ids:
            return False, "id"

        if self.config.exclude_ids and sv_id in self.config.exclude_ids:
            return False, "exclude_id"

        # SVTYPE filters
        if self.config.svtypes and svtype not in self.config.svtypes:
            return False, "svtype"

        if self.config.exclude_svtypes and svtype in self.config.exclude_svtypes:
            return False, "exclude_svtype"

        # FILTER filters
        if self.config.pass_only and filter_value != "PASS":
            return False, "pass_only"

        if self.config.filters and filter_value not in self.config.filters:
            return False, "filter"

        if self.config.exclude_filters and filter_value in self.config.exclude_filters:
            return False, "exclude_filter"

        filter_tokens = split_filter_tokens(filter_value)

        if self.config.filter_tokens and not (filter_tokens & self.config.filter_tokens):
            return False, "filter_token"

        if self.config.exclude_filter_tokens and (
            filter_tokens & self.config.exclude_filter_tokens
        ):
            return False, "exclude_filter_token"

        # Contig filters
        record_contigs = self._record_contigs(chrom, svtype, info)

        if self.config.contigs:
            if not any(contig_matches(c, self.config.contigs) for c in record_contigs):
                return False, "contig"

        if self.config.exclude_contigs:
            if any(contig_matches(c, self.config.exclude_contigs) for c in record_contigs):
                return False, "exclude_contig"

        if self.config.standard_contigs:
            if not all(is_standard_contig(c) for c in record_contigs):
                return False, "standard_contigs"

        # Size filters
        if self.config.min_size is not None or self.config.max_size is not None:
            size = svcf_event_size(svtype, pos, info)

            if size is None:
                if svtype in NON_LINEAR_SVTYPES:
                    return False, "size_not_available_for_TRA_or_BND"
                return False, "size_missing_or_invalid"

            if self.config.min_size is not None and size < self.config.min_size:
                return False, "size_below_min"

            if self.config.max_size is not None and size > self.config.max_size:
                return False, "size_above_max"

        # SUPPORT filters
        if self.config.min_support is not None or self.config.max_support is not None:
            support = safe_int(info.get("SUPPORT"))

            if support is None:
                return False, "support_missing_or_invalid"

            if self.config.min_support is not None and support < self.config.min_support:
                return False, "support_below_min"

            if self.config.max_support is not None and support > self.config.max_support:
                return False, "support_above_max"

        # Source filters
        if self._source_filters_active():
            sources, method = self._record_sources(info, sv_id)
            self.source_detection_used[method] += 1

            source_count = len(sources)

            if self.config.sources:
                if self.config.source_mode == "all":
                    if not self.config.sources.issubset(sources):
                        return False, "source"
                else:
                    if not (self.config.sources & sources):
                        return False, "source"

            if self.config.exclude_sources and (self.config.exclude_sources & sources):
                return False, "exclude_source"

            if (
                self.config.min_source_count is not None
                and source_count < self.config.min_source_count
            ):
                return False, "source_count_below_min"

            if (
                self.config.max_source_count is not None
                and source_count > self.config.max_source_count
            ):
                return False, "source_count_above_max"

        # QUAL filters
        if self.config.min_qual is not None or self.config.max_qual is not None:
            qual_value = safe_float(qual)

            if qual_value is None:
                return False, "qual_missing_or_invalid"

            if self.config.min_qual is not None and qual_value < self.config.min_qual:
                return False, "qual_below_min"

            if self.config.max_qual is not None and qual_value > self.config.max_qual:
                return False, "qual_above_max"

        # AF filters
        if self.config.min_af is not None or self.config.max_af is not None:
            af_value = safe_float(info.get("AF"))

            if af_value is None:
                return False, "af_missing_or_invalid"

            if self.config.min_af is not None and af_value < self.config.min_af:
                return False, "af_below_min"

            if self.config.max_af is not None and af_value > self.config.max_af:
                return False, "af_above_max"

        return True, "kept"

    def _record_contigs(self, chrom: str, svtype: str, info: dict) -> list[str]:
        contigs = [chrom]

        chr2 = info.get("CHR2")
        if svtype in NON_LINEAR_SVTYPES and chr2 not in (None, "", ".", True):
            contigs.append(str(chr2))

        return contigs

    def _source_filters_active(self) -> bool:
        return any(
            [
                self.config.sources,
                self.config.exclude_sources,
                self.config.min_source_count is not None,
                self.config.max_source_count is not None,
            ]
        )

    def _record_sources(self, info: dict, sv_id: str) -> tuple[set[str], str]:
        """Infer record-level caller/source support.

        Priority:
            1. INFO/SOURCES
            2. record ID prefix, e.g. pbsv.INS.123 -> pbsv
            3. single trailing column name in #CHROM header

        This is record-level provenance detection. It does not delete columns.
        """
        sources_value = info.get("SOURCES")
        if sources_value not in (None, "", ".", True):
            sources = {
                item.strip().lower()
                for item in str(sources_value).split(",")
                if item.strip()
            }
            return sources, "INFO/SOURCES"

        if "." in sv_id:
            prefix = sv_id.split(".", 1)[0].strip()
            if prefix:
                return {prefix.lower()}, "ID_prefix"

        if len(self.header_trailing_columns) == 1:
            column = self.header_trailing_columns[0].strip()
            if column:
                return {column.lower()}, "single_header_column"

        return set(), "unresolved"

    def summary(self) -> dict:
        removed = self.input_records - self.output_records

        return {
            "tool": "octopusv filter",
            "mode": "row_filter",
            "preserves_header": True,
            "preserves_columns": True,
            "input": self.config.input_file,
            "output": None if self.config.dry_run else self.config.output_file,
            "dry_run": self.config.dry_run,
            "input_records": self.input_records,
            "output_records": self.output_records,
            "removed_records": removed,
            "kept_fraction": (
                self.output_records / self.input_records if self.input_records else 0.0
            ),
            "filters": self.config.active_filters(),
            "excluded_by_reason": dict(self.excluded_by_reason),
            "source_detection_used": dict(self.source_detection_used),
        }

    @staticmethod
    def summary_to_json(summary: dict) -> str:
        return json.dumps(summary, indent=2)

    @staticmethod
    def summary_to_text(summary: dict) -> str:
        lines = [
            "SVCF filter summary",
            f"Input: {summary['input']}",
            f"Output: {summary['output']}",
            f"Dry run: {summary['dry_run']}",
            f"Input records: {summary['input_records']}",
            f"Output records: {summary['output_records']}",
            f"Removed records: {summary['removed_records']}",
            f"Kept fraction: {summary['kept_fraction']:.4f}",
            f"Preserves columns: {summary['preserves_columns']}",
        ]

        if summary["excluded_by_reason"]:
            lines.append("")
            lines.append("Excluded by reason:")
            for reason, count in sorted(summary["excluded_by_reason"].items()):
                lines.append(f"  {reason}: {count}")

        if summary["source_detection_used"]:
            lines.append("")
            lines.append("Source detection used:")
            for method, count in sorted(summary["source_detection_used"].items()):
                lines.append(f"  {method}: {count}")

        return "\n".join(lines)
