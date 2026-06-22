# src/octopusv/cli/query.py

from __future__ import annotations

from pathlib import Path
from typing import Optional

import typer

from octopusv.querying.svcf_query import (
    QueryConfig,
    SVCFQuery,
    collect_query_targets,
)


def query_svcf(
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
        help="Output matched SVCF file. Required unless --dry-run is used.",
    ),
    json_summary: bool = typer.Option(
        False,
        "--json-summary",
        help="Print an agent-friendly JSON query summary.",
    ),
    dry_run: bool = typer.Option(
        False,
        "--dry-run",
        help="Do not write an output SVCF; only report query summary.",
    ),
    region: Optional[list[str]] = typer.Option(
        None,
        "--region",
        help=(
            "Target region in 1-based closed format, e.g. chr17:7560000-7600000. "
            "Can be used multiple times."
        ),
    ),
    bed: Optional[list[Path]] = typer.Option(
        None,
        "--bed",
        exists=True,
        dir_okay=False,
        resolve_path=True,
        help=(
            "BED file containing target intervals. BED is interpreted as "
            "0-based half-open. Can be used multiple times."
        ),
    ),
    gene: Optional[list[str]] = typer.Option(
        None,
        "--gene",
        help=(
            "Gene name or gene ID to query from --gtf. Can be used multiple times. "
            "Matching is case-insensitive."
        ),
    ),
    gene_list: Optional[Path] = typer.Option(
        None,
        "--gene-list",
        exists=True,
        dir_okay=False,
        resolve_path=True,
        help="File containing gene names or gene IDs, one per line. Requires --gtf.",
    ),
    gtf: Optional[Path] = typer.Option(
        None,
        "--gtf",
        exists=True,
        dir_okay=False,
        resolve_path=True,
        help="GTF annotation file used with --gene or --gene-list.",
    ),
    gene_id_field: str = typer.Option(
        "gene_name",
        "--gene-id-field",
        help="GTF attribute used for gene matching: gene_name or gene_id.",
    ),
    flank: int = typer.Option(
        0,
        "--flank",
        help="Add this many bases to both sides of every target interval.",
    ),
    match_mode: str = typer.Option(
        "any",
        "--match-mode",
        help=(
            "Target matching mode: any, endpoint, or span. "
            "endpoint checks SV breakpoints; span checks DEL/DUP/INV intervals."
        ),
    ),
    summary_top_n: int = typer.Option(
        50,
        "--summary-top-n",
        help=(
            "Maximum number of matched records shown in the summary. "
            "Does not limit SVCF output."
        ),
    ),
    summary_target_n: int = typer.Option(
        100,
        "--summary-target-n",
        help=(
            "Maximum number of query targets shown in the summary. "
            "Does not limit target matching."
        ),
    ),
):
    """Query SVCF records by genomic targets while preserving SVCF structure."""
    if not dry_run and output_file is None:
        typer.echo(
            "Error: -o/--output-file is required unless --dry-run is used.",
            err=True,
        )
        raise typer.Exit(code=1)

    if not any([region, bed, gene, gene_list]):
        typer.echo(
            "Error: At least one target source is required: "
            "--region, --bed, --gene, or --gene-list.",
            err=True,
        )
        raise typer.Exit(code=1)

    if (gene or gene_list) and gtf is None:
        typer.echo(
            "Error: --gtf is required when using --gene or --gene-list.",
            err=True,
        )
        raise typer.Exit(code=1)

    try:
        targets, target_build_warnings = collect_query_targets(
            regions=region,
            bed_files=bed,
            genes=gene,
            gene_list_file=gene_list,
            gtf_file=gtf,
            flank=flank,
            gene_id_field=gene_id_field,
        )

        config = QueryConfig(
            input_file=str(input_file),
            output_file=str(output_file) if output_file is not None else None,
            dry_run=dry_run,
            targets=targets,
            match_mode=match_mode,
            summary_top_n=summary_top_n,
            summary_target_n=summary_target_n,
            target_build_warnings=target_build_warnings,
        )

        summary = SVCFQuery(config).run()

    except (FileNotFoundError, ValueError) as exc:
        typer.echo(f"Error: {exc}", err=True)
        raise typer.Exit(code=1) from exc

    if json_summary:
        typer.echo(SVCFQuery.summary_to_json(summary))
    else:
        typer.echo(SVCFQuery.summary_to_text(summary))


if __name__ == "__main__":
    typer.run(query_svcf)