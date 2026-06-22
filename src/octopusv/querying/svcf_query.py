"""SVCF target-based query engine.

This module implements target-based record retrieval for OctopuSV SVCF files.

Important boundary:
    query = target-based row retrieval

Therefore:
    octopusv query
        retrieves records whose breakpoints or linear SV intervals hit
        user-provided genomic targets.

It preserves all SVCF header lines and writes matched records without rewriting
INFO/FORMAT/sample/caller columns.

Coordinate convention:
    QueryTarget uses 1-based closed intervals internally.

Examples:
    VCF/SVCF region: chr1:100-200 -> start=100, end=200
    BED interval later: chr1 99 200 -> start=100, end=200

v0.1 scope:
    --region
    --flank
    --match-mode any|endpoint|span

v0.1 intentionally does not implement:
    --bed
    --gene / --gene-list / --gtf
    annotation / consequence prediction / pathogenicity ranking
"""

from __future__ import annotations

import json
import re
from collections import Counter
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

from octopusv.filtering.svcf_filter import contig_matches, parse_info, safe_int


NON_LINEAR_SVTYPES = {"TRA", "BND"}
SPAN_SVTYPES = {"DEL", "DUP", "INV"}
LEGAL_MATCH_MODES = {"any", "endpoint", "span"}


@dataclass
class QueryTarget:
    """A genomic target interval.

    Internal coordinate convention:
        1-based closed interval.

    Fields:
        target_id:
            Human-readable target label, e.g. "chr17:7560000-7600000".
        source:
            Where the target came from, e.g. "region".
        chrom:
            Target chromosome/contig.
        start:
            1-based closed start.
        end:
            1-based closed end.
        flank:
            Flank added to the original target.
        metadata:
            Optional structured metadata, such as original_start/original_end.
    """

    target_id: str
    source: str
    chrom: str
    start: int
    end: int
    flank: int = 0
    metadata: dict = field(default_factory=dict)

    def __post_init__(self) -> None:
        if self.start < 1:
            self.start = 1

        if self.end < self.start:
            raise ValueError(
                f"Invalid target interval for {self.target_id}: "
                f"{self.chrom}:{self.start}-{self.end}"
            )

    def to_dict(self) -> dict:
        """Return a JSON-friendly representation."""
        return {
            "target_id": self.target_id,
            "source": self.source,
            "chrom": self.chrom,
            "start": self.start,
            "end": self.end,
            "flank": self.flank,
            "metadata": self.metadata,
        }


@dataclass
class QueryConfig:
    """Configuration for SVCF target-based query."""

    input_file: str
    output_file: Optional[str] = None
    dry_run: bool = False
    targets: list[QueryTarget] = field(default_factory=list)
    match_mode: str = "any"
    summary_top_n: int = 50

    def __post_init__(self) -> None:
        self.match_mode = str(self.match_mode).lower()

        if self.match_mode not in LEGAL_MATCH_MODES:
            raise ValueError(
                f"Invalid match mode '{self.match_mode}'. "
                f"Allowed values: {sorted(LEGAL_MATCH_MODES)}."
            )

        if self.summary_top_n < 0:
            raise ValueError("--summary-top-n must be non-negative.")

        if not self.targets:
            raise ValueError(
                "At least one target is required. Use --region for v0.1 query."
            )


@dataclass
class Endpoint:
    """A record breakpoint endpoint."""

    chrom: str
    pos: int
    label: str  # endpoint1 or endpoint2


@dataclass
class Span:
    """A record reference span.

    v0.1 only computes spans for DEL/DUP/INV.

    INS is intentionally excluded because in OctopuSV SVCF, INS END does not
    represent a reference affected interval.
    """

    chrom: str
    start: int
    end: int


def parse_region(region: str, flank: int = 0) -> QueryTarget:
    """Parse a region string like 'chr17:7560000-7600000'.

    The input region is expected to be 1-based closed.

    Commas in coordinates are accepted:
        chr17:7,560,000-7,600,000
    """
    if flank < 0:
        raise ValueError("--flank must be non-negative.")

    text = str(region).strip().replace(",", "")
    if not text:
        raise ValueError("Empty --region value.")

    if ":" not in text or "-" not in text:
        raise ValueError(
            f"Invalid region '{region}'. Expected format: chrom:start-end."
        )

    chrom, rest = text.split(":", 1)
    start_str, end_str = rest.split("-", 1)

    chrom = chrom.strip()
    if not chrom:
        raise ValueError(f"Invalid region '{region}'. Chromosome is empty.")

    try:
        start = int(start_str)
        end = int(end_str)
    except ValueError as exc:
        raise ValueError(
            f"Invalid region '{region}'. Start/end must be integers."
        ) from exc

    if start < 1 or end < 1 or end < start:
        raise ValueError(
            f"Invalid region '{region}'. Expected 1-based closed start <= end."
        )

    original_start = start
    original_end = end

    start = max(1, start - flank)
    end = end + flank

    return QueryTarget(
        target_id=f"{chrom}:{original_start}-{original_end}",
        source="region",
        chrom=chrom,
        start=start,
        end=end,
        flank=flank,
        metadata={
            "original_start": original_start,
            "original_end": original_end,
        },
    )


def parse_regions(regions: list[str], flank: int = 0) -> list[QueryTarget]:
    """Parse multiple --region values into QueryTarget objects."""
    targets: list[QueryTarget] = []

    for region in regions:
        targets.append(parse_region(region, flank=flank))

    return targets


class SVCFQuery:
    """Retrieve SVCF records by genomic targets while preserving SVCF structure."""

    def __init__(self, config: QueryConfig):
        self.config = config

        self.input_records = 0
        self.matched_records = 0

        self.matched_by_reason = Counter()
        self.svtype_counts = Counter()
        self.skipped_or_not_matched_by_reason = Counter()

        self.target_hit_counts = Counter()
        self.top_matched_records: list[dict] = []

    def run(self) -> dict:
        """Run query and return a summary dict."""
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
                        if output_handle is not None:
                            output_handle.write(line)
                        continue

                    if not line.strip():
                        continue

                    self.input_records += 1

                    matched, record_summary = self._match_line(line)
                    if not matched:
                        continue

                    self.matched_records += 1
                    self.svtype_counts[record_summary["svtype"]] += 1

                    for reason in record_summary["match_reason"]:
                        self.matched_by_reason[reason] += 1

                    for target_id in record_summary["matched_targets"]:
                        self.target_hit_counts[target_id] += 1

                    if (
                        self.config.summary_top_n > 0
                        and len(self.top_matched_records) < self.config.summary_top_n
                    ):
                        self.top_matched_records.append(record_summary)

                    if output_handle is not None:
                        output_handle.write(line)

        finally:
            if output_handle is not None:
                output_handle.close()

        return self.summary()

    def _match_line(self, line: str) -> tuple[bool, dict]:
        """Return whether a record matches targets and its summary if matched."""
        fields = line.rstrip("\n").split("\t")
        if len(fields) < 8:
            self.skipped_or_not_matched_by_reason["malformed_record"] += 1
            return False, {}

        chrom = fields[0]
        pos = fields[1]
        sv_id = fields[2]
        alt = fields[4] if len(fields) > 4 else "."
        info = parse_info(fields[7])

        svtype = str(info.get("SVTYPE", "Unknown")).upper()

        endpoints = self._record_endpoints(
            chrom=chrom,
            pos=pos,
            alt=alt,
            svtype=svtype,
            info=info,
        )

        span, span_unavailable_reason = self._record_span(
            chrom=chrom,
            pos=pos,
            svtype=svtype,
            info=info,
        )

        endpoint_hits: dict[str, set[str]] = {}
        span_hits: dict[str, str] = {}

        if self.config.match_mode in {"any", "endpoint"}:
            endpoint_hits = self._match_endpoints(endpoints)

        if self.config.match_mode in {"any", "span"}:
            if span is None:
                # Avoid noisy counts in default "any" mode. In "span" mode,
                # users explicitly asked for span matching, so report why
                # a span was unavailable.
                if self.config.match_mode == "span" and span_unavailable_reason:
                    self.skipped_or_not_matched_by_reason[
                        span_unavailable_reason
                    ] += 1
            else:
                span_hits = self._match_span(span)

        endpoint_target_ids: set[str] = set()
        for target_ids in endpoint_hits.values():
            endpoint_target_ids.update(target_ids)

        matched_targets = sorted(endpoint_target_ids | set(span_hits.values()))

        if not matched_targets:
            return False, {}

        reasons = set()

        endpoint_reason_present = False
        span_reason_present = False

        for endpoint_label in endpoint_hits:
            reasons.add(f"{endpoint_label}_in_target")
            endpoint_reason_present = True

        if span_hits:
            reasons.add("span_overlap_target")
            span_reason_present = True

        if endpoint_reason_present and span_reason_present:
            reasons.add("endpoint_and_span")

        endpoint2 = self._get_endpoint_by_label(endpoints, "endpoint2")

        record_summary = {
            "id": sv_id,
            "svtype": svtype,
            "chrom": chrom,
            "pos": safe_int(pos),
            "chr2": endpoint2.chrom if endpoint2 is not None else None,
            "end": endpoint2.pos if endpoint2 is not None else safe_int(info.get("END")),
            "matched_targets": matched_targets,
            "match_reason": sorted(reasons),
        }

        return True, record_summary

    def _record_endpoints(
        self,
        chrom: str,
        pos: str,
        alt: str,
        svtype: str,
        info: dict,
    ) -> list[Endpoint]:
        """Return query endpoints for a record.

        endpoint1:
            Always CHROM:POS when POS is valid.

        endpoint2:
            DEL/DUP/INV:
                CHROM:END when END is valid.
            TRA/BND:
                CHR2:END when valid, otherwise BND-like ALT fallback.
            INS:
                No endpoint2 in v0.1. Treated as point-like at CHROM:POS.
        """
        endpoints: list[Endpoint] = []

        pos_i = safe_int(pos)
        if pos_i is not None:
            endpoints.append(Endpoint(chrom=chrom, pos=pos_i, label="endpoint1"))

        end_i = safe_int(info.get("END"))
        chr2 = self._clean_info_value(info.get("CHR2"))

        if svtype in NON_LINEAR_SVTYPES:
            if chr2 is not None and end_i is not None:
                endpoints.append(Endpoint(chrom=chr2, pos=end_i, label="endpoint2"))
                return endpoints

            alt_endpoint = self._endpoint_from_bnd_alt(alt)
            if alt_endpoint is not None:
                endpoints.append(
                    Endpoint(
                        chrom=alt_endpoint[0],
                        pos=alt_endpoint[1],
                        label="endpoint2",
                    )
                )
            return endpoints

        if svtype in SPAN_SVTYPES and end_i is not None:
            endpoints.append(Endpoint(chrom=chrom, pos=end_i, label="endpoint2"))

        # INS is intentionally endpoint1-only in query v0.1.
        return endpoints

    def _record_span(
        self,
        chrom: str,
        pos: str,
        svtype: str,
        info: dict,
    ) -> tuple[Optional[Span], Optional[str]]:
        """Return reference span for DEL/DUP/INV only.

        INS is intentionally excluded because INS END in OctopuSV SVCF does not
        represent a reference affected interval.
        """
        if svtype in NON_LINEAR_SVTYPES:
            return None, "span_not_available_for_TRA_or_BND"

        if svtype == "INS":
            return None, "span_not_available_for_INS"

        if svtype not in SPAN_SVTYPES:
            return None, "span_not_available_for_svtype"

        pos_i = safe_int(pos)
        if pos_i is None:
            return None, "span_missing_or_invalid"

        end_i = safe_int(info.get("END"))

        if end_i is None:
            return None, "span_missing_or_invalid"

        start = min(pos_i, end_i)
        end = max(pos_i, end_i)

        return Span(chrom=chrom, start=start, end=end), None

    def _match_endpoints(self, endpoints: list[Endpoint]) -> dict[str, set[str]]:
        """Return endpoint label -> target IDs for endpoint hits."""
        hits: dict[str, set[str]] = {}

        for endpoint in endpoints:
            for target in self.config.targets:
                if self._point_in_target(endpoint.chrom, endpoint.pos, target):
                    hits.setdefault(endpoint.label, set()).add(target.target_id)

        return hits

    def _match_span(self, span: Span) -> dict[str, str]:
        """Return target_id -> target_id for span-overlap hits."""
        hits: dict[str, str] = {}

        for target in self.config.targets:
            if not contig_matches(span.chrom, {target.chrom}):
                continue

            if self._closed_intervals_overlap(
                span.start,
                span.end,
                target.start,
                target.end,
            ):
                hits[target.target_id] = target.target_id

        return hits

    def _point_in_target(self, chrom: str, pos: int, target: QueryTarget) -> bool:
        """Return True if a point falls inside a target interval."""
        if not contig_matches(chrom, {target.chrom}):
            return False

        return target.start <= pos <= target.end

    @staticmethod
    def _closed_intervals_overlap(
        start1: int,
        end1: int,
        start2: int,
        end2: int,
    ) -> bool:
        """Return True if two 1-based closed intervals overlap."""
        return start1 <= end2 and end1 >= start2

    @staticmethod
    def _clean_info_value(value) -> Optional[str]:
        """Normalize missing INFO values."""
        if value in (None, "", ".", True):
            return None
        return str(value)

    @staticmethod
    def _endpoint_from_bnd_alt(alt: str) -> Optional[tuple[str, int]]:
        """Parse remote endpoint from a BND-like ALT field.

        Examples:
            N]chr2:12345]
            ]chr2:12345]N
            N[chr2:12345[
            [chr2:12345[N
        """
        if not alt or alt == ".":
            return None

        match = re.search(r"[\[\]]([^:\[\]]+):(\d+)[\[\]]", alt)
        if not match:
            return None

        chrom = match.group(1)
        pos = safe_int(match.group(2))
        if pos is None:
            return None

        return chrom, pos

    @staticmethod
    def _get_endpoint_by_label(
        endpoints: list[Endpoint],
        label: str,
    ) -> Optional[Endpoint]:
        for endpoint in endpoints:
            if endpoint.label == label:
                return endpoint
        return None

    def summary(self) -> dict:
        """Return an agent-friendly structured summary."""
        target_summaries = [target.to_dict() for target in self.config.targets]

        unmatched_targets = [
            target.target_id
            for target in self.config.targets
            if self.target_hit_counts[target.target_id] == 0
        ]

        return {
            "tool": "octopusv query",
            "mode": "target_query",
            "preserves_header": True,
            "preserves_columns": True,
            "input": self.config.input_file,
            "output": None if self.config.dry_run else self.config.output_file,
            "dry_run": self.config.dry_run,
            "input_records": self.input_records,
            "matched_records": self.matched_records,
            "unmatched_records": self.input_records - self.matched_records,
            "matched_fraction": (
                self.matched_records / self.input_records if self.input_records else 0.0
            ),
            "match_mode": self.config.match_mode,
            "targets": target_summaries,
            "matched_by_reason": dict(self.matched_by_reason),
            "svtype_counts": dict(self.svtype_counts),
            "top_matched_records": self.top_matched_records,
            "summary_top_n": self.config.summary_top_n,
            "unmatched_targets": unmatched_targets,
            "skipped_or_not_matched_by_reason": dict(
                self.skipped_or_not_matched_by_reason
            ),
        }

    @staticmethod
    def summary_to_json(summary: dict) -> str:
        return json.dumps(summary, indent=2)

    @staticmethod
    def summary_to_text(summary: dict) -> str:
        lines = [
            "SVCF query summary",
            f"Input: {summary['input']}",
            f"Output: {summary['output']}",
            f"Dry run: {summary['dry_run']}",
            f"Input records: {summary['input_records']}",
            f"Matched records: {summary['matched_records']}",
            f"Unmatched records: {summary['unmatched_records']}",
            f"Matched fraction: {summary['matched_fraction']:.4f}",
            f"Match mode: {summary['match_mode']}",
            f"Targets: {len(summary['targets'])}",
            f"Preserves columns: {summary['preserves_columns']}",
        ]

        if summary["matched_by_reason"]:
            lines.append("")
            lines.append("Matched by reason:")
            for reason, count in sorted(summary["matched_by_reason"].items()):
                lines.append(f"  {reason}: {count}")

        if summary["svtype_counts"]:
            lines.append("")
            lines.append("Matched SVTYPE counts:")
            for svtype, count in sorted(summary["svtype_counts"].items()):
                lines.append(f"  {svtype}: {count}")

        if summary["skipped_or_not_matched_by_reason"]:
            lines.append("")
            lines.append("Skipped or not matched by reason:")
            for reason, count in sorted(
                summary["skipped_or_not_matched_by_reason"].items()
            ):
                lines.append(f"  {reason}: {count}")

        if summary["unmatched_targets"]:
            lines.append("")
            lines.append("Unmatched targets:")
            for target_id in summary["unmatched_targets"][:20]:
                lines.append(f"  {target_id}")
            if len(summary["unmatched_targets"]) > 20:
                remaining = len(summary["unmatched_targets"]) - 20
                lines.append(f"  ... {remaining} more")

        return "\n".join(lines)