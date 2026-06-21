# src/octopusv/cli/filter.py

from __future__ import annotations

from pathlib import Path
from typing import Optional

import typer

from octopusv.filtering.svcf_filter import (
    FilterConfig,
    SVCFFilter,
    parse_csv_arg,
    read_list_file,
)


def filter_svcf(
    input_file: Path = typer.Option(
        ...,
        "--input-file",
        "-i",
        exists=True,
        dir_okay=False,
        resolve_path=True,
        help="Input SVCF file.",
    ),
    output_file: Optional[Path] = typer.Option(
        None,
        "--output-file",
        "-o",
        dir_okay=False,
        resolve_path=True,
        help="Output filtered SVCF file. Required unless --dry-run is used.",
    ),
    json_summary: bool = typer.Option(
        False,
        "--json-summary",
        help="Print an agent-friendly JSON filtering summary.",
    ),
    dry_run: bool = typer.Option(
        False,
        "--dry-run",
        help="Do not write an output SVCF; only report filtering summary.",
    ),
    svtype: Optional[str] = typer.Option(
        None,
        "--svtype",
        help="Comma-separated SVTYPE include list, e.g. DEL,DUP,INV,INS,TRA.",
    ),
    exclude_svtype: Optional[str] = typer.Option(
        None,
        "--exclude-svtype",
        help="Comma-separated SVTYPE exclusion list.",
    ),
    pass_only: bool = typer.Option(
        False,
        "--pass-only",
        help="Keep only records with FILTER exactly equal to PASS.",
    ),
    filter_value: Optional[str] = typer.Option(
        None,
        "--filter",
        help="Comma-separated exact FILTER values to keep.",
    ),
    exclude_filter: Optional[str] = typer.Option(
        None,
        "--exclude-filter",
        help="Comma-separated exact FILTER values to exclude.",
    ),
    filter_token: Optional[str] = typer.Option(
        None,
        "--filter-token",
        help="Comma-separated FILTER tokens to keep after splitting FILTER by ';'.",
    ),
    exclude_filter_token: Optional[str] = typer.Option(
        None,
        "--exclude-filter-token",
        help="Comma-separated FILTER tokens to exclude after splitting FILTER by ';'.",
    ),
    min_size: Optional[int] = typer.Option(
        None,
        "--min-size",
        help=(
            "Minimum SVCF-defined event size. DEL/DUP/INV/INS use SVLEN or END-POS; "
            "TRA/BND do not pass size filters."
        ),
    ),
    max_size: Optional[int] = typer.Option(
        None,
        "--max-size",
        help=(
            "Maximum SVCF-defined event size. DEL/DUP/INV/INS use SVLEN or END-POS; "
            "TRA/BND do not pass size filters."
        ),
    ),
    min_support: Optional[int] = typer.Option(
        None,
        "--min-support",
        help="Minimum INFO/SUPPORT read support. Missing or '.' SUPPORT does not pass.",
    ),
    max_support: Optional[int] = typer.Option(
        None,
        "--max-support",
        help="Maximum INFO/SUPPORT read support. Missing or '.' SUPPORT does not pass.",
    ),
    source: Optional[str] = typer.Option(
        None,
        "--source",
        help=(
            "Comma-separated record-level source/caller support include list. "
            "Uses INFO/SOURCES first, then ID prefix and single-header-column fallback."
        ),
    ),
    exclude_source: Optional[str] = typer.Option(
        None,
        "--exclude-source",
        help="Comma-separated record-level source/caller support exclusion list.",
    ),
    source_mode: str = typer.Option(
        "any",
        "--source-mode",
        help="How to apply --source: any or all.",
    ),
    min_source_count: Optional[int] = typer.Option(
        None,
        "--min-source-count",
        help="Keep records supported by at least this many sources/callers.",
    ),
    max_source_count: Optional[int] = typer.Option(
        None,
        "--max-source-count",
        help="Keep records supported by at most this many sources/callers.",
    ),
    record_id: Optional[str] = typer.Option(
        None,
        "--id",
        help="Comma-separated record IDs to keep.",
    ),
    id_list: Optional[Path] = typer.Option(
        None,
        "--id-list",
        exists=True,
        dir_okay=False,
        resolve_path=True,
        help="File containing record IDs to keep, one per line.",
    ),
    exclude_id: Optional[str] = typer.Option(
        None,
        "--exclude-id",
        help="Comma-separated record IDs to exclude.",
    ),
    exclude_id_list: Optional[Path] = typer.Option(
        None,
        "--exclude-id-list",
        exists=True,
        dir_okay=False,
        resolve_path=True,
        help="File containing record IDs to exclude, one per line.",
    ),
    contig: Optional[str] = typer.Option(
        None,
        "--contig",
        help=(
            "Comma-separated contigs to keep. For TRA/BND, either endpoint matching "
            "the contig is sufficient."
        ),
    ),
    exclude_contig: Optional[str] = typer.Option(
        None,
        "--exclude-contig",
        help=(
            "Comma-separated contigs to exclude. For TRA/BND, either endpoint matching "
            "causes exclusion."
        ),
    ),
    standard_contigs: bool = typer.Option(
        False,
        "--standard-contigs",
        help=(
            "Keep only records whose relevant contigs are standard chromosomes "
            "1-22/X/Y/MT or chr1-chr22/chrX/chrY/chrM. "
            "For TRA/BND, both endpoints must be standard."
        ),
    ),
    min_qual: Optional[float] = typer.Option(
        None,
        "--min-qual",
        help="Minimum numeric QUAL. Missing or '.' QUAL does not pass.",
    ),
    max_qual: Optional[float] = typer.Option(
        None,
        "--max-qual",
        help="Maximum numeric QUAL. Missing or '.' QUAL does not pass.",
    ),
    min_af: Optional[float] = typer.Option(
        None,
        "--min-af",
        help="Minimum INFO/AF. Missing or '.' AF does not pass.",
    ),
    max_af: Optional[float] = typer.Option(
        None,
        "--max-af",
        help="Maximum INFO/AF. Missing or '.' AF does not pass.",
    ),
):
    """Filter SVCF records by SV-level attributes while preserving SVCF structure."""
    if not dry_run and output_file is None:
        typer.echo(
            "Error: -o/--output-file is required unless --dry-run is used.",
            err=True,
        )
        raise typer.Exit(code=1)

    ids = parse_csv_arg(record_id) | read_list_file(id_list)
    exclude_ids = parse_csv_arg(exclude_id) | read_list_file(exclude_id_list)

    config = FilterConfig(
        input_file=str(input_file),
        output_file=str(output_file) if output_file is not None else None,
        dry_run=dry_run,
        svtypes=parse_csv_arg(svtype),
        exclude_svtypes=parse_csv_arg(exclude_svtype),
        pass_only=pass_only,
        filters=parse_csv_arg(filter_value),
        exclude_filters=parse_csv_arg(exclude_filter),
        filter_tokens=parse_csv_arg(filter_token),
        exclude_filter_tokens=parse_csv_arg(exclude_filter_token),
        min_size=min_size,
        max_size=max_size,
        min_support=min_support,
        max_support=max_support,
        sources=parse_csv_arg(source),
        exclude_sources=parse_csv_arg(exclude_source),
        source_mode=source_mode,
        min_source_count=min_source_count,
        max_source_count=max_source_count,
        ids=ids,
        exclude_ids=exclude_ids,
        contigs=parse_csv_arg(contig),
        exclude_contigs=parse_csv_arg(exclude_contig),
        standard_contigs=standard_contigs,
        min_qual=min_qual,
        max_qual=max_qual,
        min_af=min_af,
        max_af=max_af,
    )

    try:
        summary = SVCFFilter(config).run()
    except (FileNotFoundError, ValueError) as exc:
        typer.echo(f"Error: {exc}", err=True)
        raise typer.Exit(code=1) from exc

    if json_summary:
        typer.echo(SVCFFilter.summary_to_json(summary))
    else:
        typer.echo(SVCFFilter.summary_to_text(summary))


if __name__ == "__main__":
    typer.run(filter_svcf)


