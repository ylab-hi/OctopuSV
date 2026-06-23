from __future__ import annotations

import json
import re
from collections import Counter
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional


# Standard BND ALT pattern fragment. We only rewrite the mate contig and keep
# the original brackets/sequence untouched.
BND_MATE_PATTERN = re.compile(r"([\[\]])([^:\[\]]+):(\d+)([\[\]])")

# CO in OctopuSV is usually "chrom_pos-chrom_pos" and is the final FORMAT
# field. We only rewrite CO values that cleanly use standard contigs.
CO_STANDARD_PATTERN = re.compile(
    r"^(?P<c1>(?:chr)?(?:[1-9]|1[0-9]|2[0-2]|X|Y|M|MT))_"
    r"(?P<p1>\d+)-"
    r"(?P<c2>(?:chr)?(?:[1-9]|1[0-9]|2[0-2]|X|Y|M|MT))_"
    r"(?P<p2>\d+)$",
    re.IGNORECASE,
)

CONTIG_HEADER_ID_PATTERN = re.compile(r"(##contig=<[^>]*\bID=)([^,>]+)(.*>)")


@dataclass
class ContigNormalizationSummary:
    input_file: str
    output_file: Optional[str]
    style: str
    dry_run: bool = False
    records_processed: int = 0
    records_written: int = 0
    header_contigs_seen: int = 0
    header_contigs_updated: int = 0
    chrom_updated: int = 0
    chr2_updated: int = 0
    bnd_alt_updated: int = 0
    co_updated: int = 0
    co_unmodified_unrecognized_format: int = 0
    malformed_records: int = 0
    standard_mappings: Counter = field(default_factory=Counter)
    untouched_nonstandard_contigs: Counter = field(default_factory=Counter)
    warnings: list[str] = field(default_factory=list)

    def add_mapping(self, old: str, new: str) -> None:
        if old != new:
            self.standard_mappings[f"{old}->{new}"] += 1

    def add_nonstandard(self, contig: str) -> None:
        if contig not in (None, "", "."):
            self.untouched_nonstandard_contigs[str(contig)] += 1

    def to_dict(self) -> dict:
        return {
            "input_file": self.input_file,
            "output_file": self.output_file,
            "style": self.style,
            "dry_run": self.dry_run,
            "records_processed": self.records_processed,
            "records_written": self.records_written,
            "fields_updated": {
                "header_contigs": self.header_contigs_updated,
                "CHROM": self.chrom_updated,
                "INFO_CHR2": self.chr2_updated,
                "BND_ALT": self.bnd_alt_updated,
                "FORMAT_CO": self.co_updated,
            },
            "header_contigs_seen": self.header_contigs_seen,
            "co_unmodified_unrecognized_format": self.co_unmodified_unrecognized_format,
            "malformed_records": self.malformed_records,
            "standard_mappings": dict(sorted(self.standard_mappings.items())),
            "untouched_nonstandard_contigs": dict(
                sorted(self.untouched_nonstandard_contigs.items())
            ),
            "warnings": self.warnings,
            "note": (
                "This command normalizes standard contig names only. "
                "It does not perform genome-build liftover or coordinate conversion."
            ),
        }

    def to_json(self) -> str:
        return json.dumps(self.to_dict(), indent=2, sort_keys=True)

    def to_text(self) -> str:
        lines = [
            "SVCF contig normalization summary",
            f"Input: {self.input_file}",
            f"Output: {self.output_file if self.output_file else '.'}",
            f"Style: {self.style}",
            f"Dry run: {self.dry_run}",
            "",
            f"Records processed: {self.records_processed}",
            f"Records written: {self.records_written}",
            "",
            "Fields updated:",
            f"  header ##contig IDs: {self.header_contigs_updated}",
            f"  CHROM: {self.chrom_updated}",
            f"  INFO/CHR2: {self.chr2_updated}",
            f"  BND ALT mate contigs: {self.bnd_alt_updated}",
            f"  FORMAT/CO: {self.co_updated}",
        ]

        if self.co_unmodified_unrecognized_format:
            lines.append(
                "  FORMAT/CO left unchanged due to unrecognized format: "
                f"{self.co_unmodified_unrecognized_format}"
            )
        if self.malformed_records:
            lines.append(f"Malformed records skipped from normalization: {self.malformed_records}")

        lines.extend(["", "Standard contig mappings applied:"])
        if self.standard_mappings:
            for mapping, count in sorted(self.standard_mappings.items()):
                lines.append(f"  {mapping}: {count}")
        else:
            lines.append("  none")

        lines.extend(["", "Non-standard contigs left unchanged:"])
        if self.untouched_nonstandard_contigs:
            for contig, count in sorted(self.untouched_nonstandard_contigs.items()):
                lines.append(f"  {contig}: {count}")
        else:
            lines.append("  none")

        if self.warnings:
            lines.extend(["", "Warnings:"])
            for warning in self.warnings:
                lines.append(f"  - {warning}")

        lines.extend(
            [
                "",
                "Note: This command normalizes standard contig names only.",
                "It does not perform genome-build liftover or coordinate conversion.",
            ]
        )
        return "\n".join(lines)


def standard_contig_key(contig: str | None) -> Optional[str]:
    """Return canonical key 1-22/X/Y/MT for standard human contigs.

    Non-standard contigs return None and must be preserved unchanged.
    """
    if contig is None:
        return None

    c = str(contig).strip()
    if not c or c == ".":
        return None

    if c.lower().startswith("chr"):
        c = c[3:]

    c_upper = c.upper()
    if c_upper in {"M", "MT"}:
        return "MT"
    if c_upper in {"X", "Y"}:
        return c_upper

    if c.isdigit():
        n = int(c)
        if 1 <= n <= 22:
            return str(n)

    return None


def contig_to_style(contig: str | None, style: str) -> tuple[str | None, bool]:
    """Return (new_contig, is_standard).

    Non-standard contigs are returned unchanged with is_standard=False.
    """
    if contig is None:
        return None, False

    key = standard_contig_key(contig)
    if key is None:
        return str(contig), False

    if style == "chr":
        return "chrM" if key == "MT" else f"chr{key}", True
    if style == "nochr":
        return "MT" if key == "MT" else key, True

    raise ValueError("style must be 'chr' or 'nochr'")


def parse_info_items(info: str) -> list[tuple[str, Optional[str]]]:
    if info in (None, "", "."):
        return []
    items: list[tuple[str, Optional[str]]] = []
    for raw_item in str(info).split(";"):
        if raw_item == "":
            continue
        if "=" in raw_item:
            key, value = raw_item.split("=", 1)
            items.append((key, value))
        else:
            items.append((raw_item, None))
    return items


def format_info_items(items: list[tuple[str, Optional[str]]]) -> str:
    if not items:
        return "."
    out = []
    for key, value in items:
        if value is None:
            out.append(key)
        else:
            out.append(f"{key}={value}")
    return ";".join(out)


def normalize_info_chr2(
    info: str,
    style: str,
    summary: ContigNormalizationSummary,
) -> str:
    items = parse_info_items(info)
    changed = False
    new_items: list[tuple[str, Optional[str]]] = []

    for key, value in items:
        if key == "CHR2" and value not in (None, "", "."):
            new_value, is_standard = contig_to_style(value, style)
            if is_standard:
                if new_value != value:
                    summary.chr2_updated += 1
                    summary.add_mapping(value, new_value)
                    changed = True
                new_items.append((key, new_value))
            else:
                summary.add_nonstandard(value)
                new_items.append((key, value))
        else:
            new_items.append((key, value))

    return format_info_items(new_items) if changed else info


def normalize_bnd_alt(
    alt: str,
    style: str,
    summary: ContigNormalizationSummary,
) -> str:
    if not alt or alt == "." or ("[" not in alt and "]" not in alt):
        return alt

    changed_any = False

    def repl(match: re.Match) -> str:
        nonlocal changed_any
        left_bracket, contig, pos, right_bracket = match.groups()
        new_contig, is_standard = contig_to_style(contig, style)
        if is_standard:
            if new_contig != contig:
                changed_any = True
                summary.add_mapping(contig, new_contig)
            return f"{left_bracket}{new_contig}:{pos}{right_bracket}"

        summary.add_nonstandard(contig)
        return match.group(0)

    new_alt = BND_MATE_PATTERN.sub(repl, alt)
    if changed_any:
        summary.bnd_alt_updated += 1
    return new_alt


def normalize_co_value(
    co: str,
    style: str,
    summary: ContigNormalizationSummary,
) -> str:
    if co in (None, "", "."):
        return co

    match = CO_STANDARD_PATTERN.fullmatch(str(co))
    if not match:
        # CO may be non-standard, missing, or a future format. Keep it unchanged.
        if "_" in str(co) and "-" in str(co):
            summary.co_unmodified_unrecognized_format += 1
        return co

    c1 = match.group("c1")
    p1 = match.group("p1")
    c2 = match.group("c2")
    p2 = match.group("p2")

    new_c1, is_standard_1 = contig_to_style(c1, style)
    new_c2, is_standard_2 = contig_to_style(c2, style)

    if not (is_standard_1 and is_standard_2):
        if not is_standard_1:
            summary.add_nonstandard(c1)
        if not is_standard_2:
            summary.add_nonstandard(c2)
        return co

    if new_c1 != c1:
        summary.add_mapping(c1, new_c1)
    if new_c2 != c2:
        summary.add_mapping(c2, new_c2)

    new_co = f"{new_c1}_{p1}-{new_c2}_{p2}"
    if new_co != co:
        summary.co_updated += 1
    return new_co


def normalize_sample_co(
    sample_col: str,
    format_col: str,
    style: str,
    summary: ContigNormalizationSummary,
) -> str:
    keys = format_col.split(":") if format_col else []
    if "CO" not in keys:
        return sample_col

    co_idx = keys.index("CO")

    # OctopuSV writes CO as the final FORMAT field. This right-split avoids
    # accidentally splitting BND ALT strings that contain ':' earlier in the sample.
    if co_idx == len(keys) - 1:
        if ":" not in sample_col:
            return sample_col
        prefix, co = sample_col.rsplit(":", 1)
        new_co = normalize_co_value(co, style, summary)
        return f"{prefix}:{new_co}" if new_co != co else sample_col

    # Conservative fallback: only use normal split if the field count matches.
    parts = sample_col.split(":")
    if len(parts) != len(keys):
        summary.warnings.append(
            "FORMAT/CO was not the final FORMAT field in at least one record; "
            "some CO values were left unchanged."
        )
        return sample_col

    old_co = parts[co_idx]
    parts[co_idx] = normalize_co_value(old_co, style, summary)
    return ":".join(parts)


def normalize_contig_header_line(
    line: str,
    style: str,
    summary: ContigNormalizationSummary,
) -> str:
    match = CONTIG_HEADER_ID_PATTERN.match(line.rstrip("\n"))
    if not match:
        return line

    summary.header_contigs_seen += 1
    prefix, contig, suffix = match.groups()
    new_contig, is_standard = contig_to_style(contig, style)

    if not is_standard:
        summary.add_nonstandard(contig)
        return line

    if new_contig != contig:
        summary.header_contigs_updated += 1
        summary.add_mapping(contig, new_contig)
        new_line = f"{prefix}{new_contig}{suffix}"
        return new_line + "\n" if line.endswith("\n") else new_line

    return line


def is_stale_audit_header(line: str) -> bool:
    return line.startswith("##OctopuSV_contig_normalization")


def audit_headers(style: str) -> list[str]:
    return [
        f"##OctopuSV_contig_normalization_style={style}\n",
        "##OctopuSV_contig_normalization_scope=standard_contigs_only\n",
        (
            "##OctopuSV_contig_normalization_note="
            "contig names normalized only; no genome-build liftover or coordinate conversion performed\n"
        ),
    ]


def normalize_record_line(
    line: str,
    style: str,
    summary: ContigNormalizationSummary,
) -> str:
    raw = line.rstrip("\n")
    cols = raw.split("\t")

    if len(cols) < 8:
        summary.malformed_records += 1
        return line

    summary.records_processed += 1

    # CHROM
    old_chrom = cols[0]
    new_chrom, is_standard = contig_to_style(old_chrom, style)
    if is_standard:
        if new_chrom != old_chrom:
            cols[0] = new_chrom
            summary.chrom_updated += 1
            summary.add_mapping(old_chrom, new_chrom)
    else:
        summary.add_nonstandard(old_chrom)

    # ALT BND mate contig, if present.
    if len(cols) >= 5:
        cols[4] = normalize_bnd_alt(cols[4], style, summary)

    # INFO/CHR2
    cols[7] = normalize_info_chr2(cols[7], style, summary)

    # FORMAT/CO across all sample/evidence columns.
    if len(cols) >= 10:
        format_col = cols[8]
        for i in range(9, len(cols)):
            cols[i] = normalize_sample_co(cols[i], format_col, style, summary)

    return "\t".join(cols) + "\n"


class SVCFContigNormalizer:
    """Normalize standard contig names in an SVCF file without changing coordinates."""

    def __init__(
        self,
        input_file: str | Path,
        output_file: str | Path | None,
        style: str,
        dry_run: bool = False,
    ):
        if style not in {"chr", "nochr"}:
            raise ValueError("style must be 'chr' or 'nochr'")

        self.input_file = Path(input_file)
        self.output_file = Path(output_file) if output_file is not None else None
        self.style = style
        self.dry_run = dry_run

    def run(self) -> ContigNormalizationSummary:
        if not self.input_file.exists():
            raise FileNotFoundError(f"Input file not found: {self.input_file}")
        if not self.dry_run and self.output_file is None:
            raise ValueError("output_file is required unless dry_run=True")

        summary = ContigNormalizationSummary(
            input_file=str(self.input_file),
            output_file=str(self.output_file) if self.output_file else None,
            style=self.style,
            dry_run=self.dry_run,
        )

        out_handle = None
        try:
            if not self.dry_run:
                assert self.output_file is not None
                self.output_file.parent.mkdir(parents=True, exist_ok=True)
                out_handle = self.output_file.open("w")

            with self.input_file.open() as in_handle:
                audit_written = False

                for line in in_handle:
                    if line.startswith("##"):
                        if is_stale_audit_header(line):
                            continue
                        new_line = normalize_contig_header_line(line, self.style, summary)
                        if out_handle is not None:
                            out_handle.write(new_line)
                        continue

                    if line.startswith("#CHROM"):
                        if not audit_written and out_handle is not None:
                            for audit_line in audit_headers(self.style):
                                out_handle.write(audit_line)
                        audit_written = True
                        if out_handle is not None:
                            out_handle.write(line)
                        continue

                    new_line = normalize_record_line(line, self.style, summary)
                    summary.records_written += 1
                    if out_handle is not None:
                        out_handle.write(new_line)

        finally:
            if out_handle is not None:
                out_handle.close()

        # Deduplicate repeated warning text while preserving order.
        if summary.warnings:
            seen = set()
            deduped = []
            for warning in summary.warnings:
                if warning not in seen:
                    seen.add(warning)
                    deduped.append(warning)
            summary.warnings = deduped

        return summary