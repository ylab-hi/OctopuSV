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

Coordinate inputs:
    --region:
        1-based closed, e.g. chr1:100-200 -> start=100, end=200

    --bed:
        0-based half-open, e.g. chr1 99 200 -> start=100, end=200

    --gtf:
        1-based closed already. Do not shift GTF start/end.
"""

from __future__ import annotations

import json
import re
from collections import Counter, defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

from octopusv.filtering.svcf_filter import (
    normalize_standard_contig,
    parse_info,
    safe_int,
)


NON_LINEAR_SVTYPES = {"TRA", "BND"}
SPAN_SVTYPES = {"DEL", "DUP", "INV"}
LEGAL_MATCH_MODES = {"any", "endpoint", "span"}
LEGAL_GENE_ID_FIELDS = {"gene_name", "gene_id"}


@dataclass
class QueryTarget:
    """A genomic target interval.

    Internal coordinate convention:
        1-based closed interval.
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
    summary_target_n: int = 100
    target_build_warnings: list[str] = field(default_factory=list)

    def __post_init__(self) -> None:
        self.match_mode = str(self.match_mode).lower()

        if self.match_mode not in LEGAL_MATCH_MODES:
            raise ValueError(
                f"Invalid match mode '{self.match_mode}'. "
                f"Allowed values: {sorted(LEGAL_MATCH_MODES)}."
            )

        if self.summary_top_n < 0:
            raise ValueError("--summary-top-n must be non-negative.")

        if self.summary_target_n < 0:
            raise ValueError("--summary-target-n must be non-negative.")

        if not self.targets:
            raise ValueError("No query targets were resolved.")


@dataclass
class Endpoint:
    """A record breakpoint endpoint."""

    chrom: str
    pos: int
    label: str  # endpoint1 or endpoint2


@dataclass
class Span:
    """A record reference span.

    v0.1 computes spans only for DEL/DUP/INV.

    INS is intentionally excluded because OctopuSV INS END does not represent
    a reference affected interval.
    """

    chrom: str
    start: int
    end: int


def _contig_key(contig: str) -> str:
    """Return a normalized contig key for target bucketing.

    Standard human contigs collapse chr-prefixed and non-prefixed aliases:
        chr7, 7 -> STD:7
        chrM, MT -> STD:MT

    Nonstandard contigs use case-insensitive exact matching:
        GL000220.1 -> RAW:gl000220.1
    """
    norm = normalize_standard_contig(contig)
    if norm is not None:
        return f"STD:{norm}"
    return f"RAW:{str(contig).lower()}"


def _apply_flank(start: int, end: int, flank: int) -> tuple[int, int]:
    if flank < 0:
        raise ValueError("--flank must be non-negative.")
    return max(1, start - flank), end + flank


def _make_unique_target_ids(targets: list[QueryTarget]) -> list[QueryTarget]:
    """Ensure target_id values are unique.

    Summary bookkeeping uses target_id as a key. Duplicate BED names or repeated
    gene records are therefore made unique while keeping the original label in
    metadata.
    """
    seen = Counter()

    for target in targets:
        base_id = target.target_id
        seen[base_id] += 1

        if seen[base_id] > 1:
            target.metadata["original_target_id"] = base_id
            target.target_id = f"{base_id}#{seen[base_id]}"

    return targets


def parse_region(region: str, flank: int = 0) -> QueryTarget:
    """Parse a region string like 'chr17:7560000-7600000'.

    The input region is expected to be 1-based closed.
    Commas in coordinates are accepted.
    """
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
        original_start = int(start_str)
        original_end = int(end_str)
    except ValueError as exc:
        raise ValueError(
            f"Invalid region '{region}'. Start/end must be integers."
        ) from exc

    if original_start < 1 or original_end < 1 or original_end < original_start:
        raise ValueError(
            f"Invalid region '{region}'. Expected 1-based closed start <= end."
        )

    start, end = _apply_flank(original_start, original_end, flank)

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


def parse_regions(regions: Optional[list[str]], flank: int = 0) -> list[QueryTarget]:
    """Parse multiple --region values into QueryTarget objects."""
    if not regions:
        return []

    return [parse_region(region, flank=flank) for region in regions]


def parse_bed_file(path: Path | str, flank: int = 0) -> list[QueryTarget]:
    """Parse BED targets.

    BED is 0-based half-open. Internal QueryTarget is 1-based closed:
        BED start0 -> start1 = start0 + 1
        BED end0   -> end1   = end0
    """
    bed_path = Path(path)
    targets: list[QueryTarget] = []

    with bed_path.open() as handle:
        for line_number, line in enumerate(handle, start=1):
            raw = line.rstrip("\n")
            stripped = raw.strip()

            if not stripped or stripped.startswith("#"):
                continue
            if stripped.startswith("track ") or stripped.startswith("browser "):
                continue

            fields = raw.split("\t")
            if len(fields) < 3:
                raise ValueError(
                    f"Invalid BED line {line_number} in {bed_path}: expected >=3 columns."
                )

            chrom = fields[0].strip()
            if not chrom:
                raise ValueError(
                    f"Invalid BED line {line_number} in {bed_path}: empty chrom."
                )

            try:
                bed_start = int(fields[1])
                bed_end = int(fields[2])
            except ValueError as exc:
                raise ValueError(
                    f"Invalid BED line {line_number} in {bed_path}: "
                    "start/end must be integers."
                ) from exc

            if bed_start < 0 or bed_end <= bed_start:
                raise ValueError(
                    f"Invalid BED line {line_number} in {bed_path}: "
                    "expected 0-based half-open start < end."
                )

            original_start = bed_start + 1
            original_end = bed_end
            start, end = _apply_flank(original_start, original_end, flank)

            name = fields[3].strip() if len(fields) >= 4 and fields[3].strip() else None
            target_id = name if name and name != "." else f"{chrom}:{original_start}-{original_end}"

            targets.append(
                QueryTarget(
                    target_id=target_id,
                    source="bed",
                    chrom=chrom,
                    start=start,
                    end=end,
                    flank=flank,
                    metadata={
                        "bed_file": str(bed_path),
                        "line_number": line_number,
                        "bed_start_0based": bed_start,
                        "bed_end_0based": bed_end,
                        "original_start": original_start,
                        "original_end": original_end,
                        "name": name,
                    },
                )
            )

    return targets


def read_gene_list_file(path: Optional[Path | str]) -> list[str]:
    """Read one gene name/ID per line, ignoring blanks and comments."""
    if path is None:
        return []

    genes: list[str] = []
    with Path(path).open() as handle:
        for line in handle:
            value = line.strip()
            if not value or value.startswith("#"):
                continue
            genes.append(value)

    return genes


def parse_gtf_attributes(attribute_string: str) -> dict[str, str]:
    """Parse GTF/GFF-like attributes into a dict.

    Supports common GTF syntax:
        gene_id "ENSG..."; gene_name "TP53";

    Also tolerates key=value items.
    """
    attrs: dict[str, str] = {}

    for item in attribute_string.strip().strip(";").split(";"):
        item = item.strip()
        if not item:
            continue

        if "=" in item and " " not in item.split("=", 1)[0]:
            key, value = item.split("=", 1)
            attrs[key.strip()] = value.strip().strip('"')
            continue

        if " " in item:
            key, value = item.split(" ", 1)
            attrs[key.strip()] = value.strip().strip('"')

    return attrs


def parse_gtf_gene_targets(
    gtf_file: Path | str,
    genes: list[str],
    gene_id_field: str = "gene_name",
    flank: int = 0,
) -> tuple[list[QueryTarget], list[str]]:
    """Parse gene-body targets from a GTF file.

    GTF coordinates are already 1-based closed, so start/end are not shifted.
    Gene matching is case-insensitive by default.
    """
    if gene_id_field not in LEGAL_GENE_ID_FIELDS:
        raise ValueError(
            f"Invalid --gene-id-field '{gene_id_field}'. "
            f"Allowed values: {sorted(LEGAL_GENE_ID_FIELDS)}."
        )

    requested = [g.strip() for g in genes if g and g.strip()]
    requested_norm_to_original = {g.lower(): g for g in requested}

    if not requested_norm_to_original:
        return [], []

    gtf_path = Path(gtf_file)
    targets: list[QueryTarget] = []
    found_norm: set[str] = set()

    with gtf_path.open() as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line.strip() or line.startswith("#"):
                continue

            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9:
                continue

            chrom = fields[0]
            feature = fields[2]
            if feature != "gene":
                continue

            try:
                start = int(fields[3])
                end = int(fields[4])
            except ValueError:
                continue

            if start < 1 or end < start:
                continue

            attrs = parse_gtf_attributes(fields[8])
            gene_value = attrs.get(gene_id_field)
            if not gene_value:
                continue

            gene_norm = gene_value.lower()
            if gene_norm not in requested_norm_to_original:
                continue

            found_norm.add(gene_norm)

            query_start, query_end = _apply_flank(start, end, flank)

            gene_name = attrs.get("gene_name")
            gene_id = attrs.get("gene_id")

            targets.append(
                QueryTarget(
                    target_id=gene_value,
                    source="gtf",
                    chrom=chrom,
                    start=query_start,
                    end=query_end,
                    flank=flank,
                    metadata={
                        "gtf_file": str(gtf_path),
                        "line_number": line_number,
                        "gene_id_field": gene_id_field,
                        "gene_name": gene_name,
                        "gene_id": gene_id,
                        "strand": fields[6],
                        "original_start": start,
                        "original_end": end,
                    },
                )
            )

    warnings: list[str] = []
    for gene_norm, original in requested_norm_to_original.items():
        if gene_norm not in found_norm:
            warnings.append(f"gene_not_found:{original}")

    return targets, warnings


def collect_query_targets(
    regions: Optional[list[str]] = None,
    bed_files: Optional[list[Path | str]] = None,
    genes: Optional[list[str]] = None,
    gene_list_file: Optional[Path | str] = None,
    gtf_file: Optional[Path | str] = None,
    flank: int = 0,
    gene_id_field: str = "gene_name",
) -> tuple[list[QueryTarget], list[str]]:
    """Collect all target sources into one QueryTarget list."""
    targets: list[QueryTarget] = []
    warnings: list[str] = []

    targets.extend(parse_regions(regions, flank=flank))

    if bed_files:
        for bed_file in bed_files:
            targets.extend(parse_bed_file(bed_file, flank=flank))

    all_genes: list[str] = []
    if genes:
        all_genes.extend(genes)
    all_genes.extend(read_gene_list_file(gene_list_file))

    if all_genes:
        if gtf_file is None:
            raise ValueError("--gtf is required when using --gene or --gene-list.")

        gene_targets, gene_warnings = parse_gtf_gene_targets(
            gtf_file=gtf_file,
            genes=all_genes,
            gene_id_field=gene_id_field,
            flank=flank,
        )
        targets.extend(gene_targets)
        warnings.extend(gene_warnings)

    targets = _make_unique_target_ids(targets)

    return targets, warnings


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
        self.target_source_counts = Counter(target.source for target in self.config.targets)

        self.top_matched_records: list[dict] = []
        self.targets_by_contig_key = self._build_target_buckets(self.config.targets)

    @staticmethod
    def _build_target_buckets(
        targets: list[QueryTarget],
    ) -> dict[str, list[QueryTarget]]:
        buckets: dict[str, list[QueryTarget]] = defaultdict(list)

        for target in targets:
            buckets[_contig_key(target.chrom)].append(target)

        return buckets

    def _targets_for_contig(self, contig: str) -> list[QueryTarget]:
        return self.targets_by_contig_key.get(_contig_key(contig), [])

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

        summary_chr2 = endpoint2.chrom if endpoint2 is not None else None
        summary_end = endpoint2.pos if endpoint2 is not None else safe_int(info.get("END"))

        if svtype == "INS":
            summary_chr2 = None
            summary_end = None

        record_summary = {
            "id": sv_id,
            "svtype": svtype,
            "chrom": chrom,
            "pos": safe_int(pos),
            "chr2": summary_chr2,
            "end": summary_end,
            "matched_targets": matched_targets,
            "matched_target_count": len(matched_targets),
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

        return endpoints

    def _record_span(
        self,
        chrom: str,
        pos: str,
        svtype: str,
        info: dict,
    ) -> tuple[Optional[Span], Optional[str]]:
        """Return reference span for DEL/DUP/INV only.

        INS is intentionally excluded because OctopuSV INS END does not
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
            svlen = safe_int(info.get("SVLEN"))
            if svlen is not None:
                end_i = pos_i + abs(svlen)

        if end_i is None:
            return None, "span_missing_or_invalid"

        start = min(pos_i, end_i)
        end = max(pos_i, end_i)

        return Span(chrom=chrom, start=start, end=end), None

    def _match_endpoints(self, endpoints: list[Endpoint]) -> dict[str, set[str]]:
        """Return endpoint label -> target IDs for endpoint hits."""
        hits: dict[str, set[str]] = {}

        for endpoint in endpoints:
            for target in self._targets_for_contig(endpoint.chrom):
                if self._point_in_target(endpoint.pos, target):
                    hits.setdefault(endpoint.label, set()).add(target.target_id)

        return hits

    def _match_span(self, span: Span) -> dict[str, str]:
        """Return target_id -> target_id for span-overlap hits."""
        hits: dict[str, str] = {}

        for target in self._targets_for_contig(span.chrom):
            if self._closed_intervals_overlap(
                span.start,
                span.end,
                target.start,
                target.end,
            ):
                hits[target.target_id] = target.target_id

        return hits

    @staticmethod
    def _point_in_target(pos: int, target: QueryTarget) -> bool:
        """Return True if a point falls inside a target interval."""
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
        targets_total = len(self.config.targets)
        targets_shown = min(targets_total, self.config.summary_target_n)

        target_summaries = [
            target.to_dict()
            for target in self.config.targets[: self.config.summary_target_n]
        ]

        unmatched_all = [
            target.target_id
            for target in self.config.targets
            if self.target_hit_counts[target.target_id] == 0
        ]
        unmatched_shown = min(len(unmatched_all), self.config.summary_target_n)

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
            "match_reason_counts_are_nonexclusive": True,
            "targets_total": targets_total,
            "targets_shown": targets_shown,
            "targets_truncated": targets_shown < targets_total,
            "targets": target_summaries,
            "target_source_counts": dict(self.target_source_counts),
            "matched_by_reason": dict(self.matched_by_reason),
            "svtype_counts": dict(self.svtype_counts),
            "top_matched_records": self.top_matched_records,
            "summary_top_n": self.config.summary_top_n,
            "summary_target_n": self.config.summary_target_n,
            "unmatched_targets_total": len(unmatched_all),
            "unmatched_targets_shown": unmatched_shown,
            "unmatched_targets_truncated": unmatched_shown < len(unmatched_all),
            "unmatched_targets": unmatched_all[: self.config.summary_target_n],
            "target_build_warnings": self.config.target_build_warnings,
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
            f"Matched fraction: {summary['matched_fraction']:.6g}",
            f"Match mode: {summary['match_mode']}",
            f"Targets total: {summary['targets_total']}",
            f"Targets shown in summary: {summary['targets_shown']}",
            f"Preserves columns: {summary['preserves_columns']}",
        ]

        if summary["target_source_counts"]:
            lines.append("")
            lines.append("Target sources:")
            for source, count in sorted(summary["target_source_counts"].items()):
                lines.append(f"  {source}: {count}")

        if summary["matched_by_reason"]:
            lines.append("")
            lines.append("Matched by reason (non-exclusive):")
            for reason, count in sorted(summary["matched_by_reason"].items()):
                lines.append(f"  {reason}: {count}")

        if summary["svtype_counts"]:
            lines.append("")
            lines.append("Matched SVTYPE counts:")
            for svtype, count in sorted(summary["svtype_counts"].items()):
                lines.append(f"  {svtype}: {count}")

        if summary["target_build_warnings"]:
            lines.append("")
            lines.append("Target build warnings:")
            for warning in summary["target_build_warnings"][:20]:
                lines.append(f"  {warning}")
            if len(summary["target_build_warnings"]) > 20:
                remaining = len(summary["target_build_warnings"]) - 20
                lines.append(f"  ... {remaining} more")

        if summary["skipped_or_not_matched_by_reason"]:
            lines.append("")
            lines.append("Skipped or not matched by reason:")
            for reason, count in sorted(
                summary["skipped_or_not_matched_by_reason"].items()
            ):
                lines.append(f"  {reason}: {count}")

        if summary["unmatched_targets"]:
            lines.append("")
            lines.append(
                f"Unmatched targets "
                f"({summary['unmatched_targets_shown']}/"
                f"{summary['unmatched_targets_total']} shown):"
            )
            for target_id in summary["unmatched_targets"]:
                lines.append(f"  {target_id}")
            if summary["unmatched_targets_truncated"]:
                remaining = (
                    summary["unmatched_targets_total"]
                    - summary["unmatched_targets_shown"]
                )
                lines.append(f"  ... {remaining} more")

        return "\n".join(lines)