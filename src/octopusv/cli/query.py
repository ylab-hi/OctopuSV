# src/octopusv/cli/query.py

from __future__ import annotations

from pathlib import Path
from typing import Optional

import typer

from octopusv.querying.svcf_query import QueryConfig, SVCFQuery, parse_regions


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
):
    """Query SVCF records by genomic targets while preserving SVCF structure."""
    if not dry_run and output_file is None:
        typer.echo(
            "Error: -o/--output-file is required unless --dry-run is used.",
            err=True,
        )
        raise typer.Exit(code=1)

    if not region:
        typer.echo(
            "Error: At least one target is required. For v0.1, use --region.",
            err=True,
        )
        raise typer.Exit(code=1)

    try:
        targets = parse_regions(region, flank=flank)

        config = QueryConfig(
            input_file=str(input_file),
            output_file=str(output_file) if output_file is not None else None,
            dry_run=dry_run,
            targets=targets,
            match_mode=match_mode,
            summary_top_n=summary_top_n,
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