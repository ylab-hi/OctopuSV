from __future__ import annotations

from pathlib import Path
from typing import Optional

import typer

from octopusv.normalization.contig_normalizer import SVCFContigNormalizer


def normalize_contigs(
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
        help="Output SVCF file. Required unless --dry-run is used.",
    ),
    style: str = typer.Option(
        ...,
        "--style",
        help="Target standard-contig naming style: chr or nochr.",
    ),
    dry_run: bool = typer.Option(
        False,
        "--dry-run",
        help="Do not write an output SVCF; only report what would be normalized.",
    ),
    json_summary: bool = typer.Option(
        False,
        "--json-summary",
        help="Print an agent-friendly JSON normalization summary.",
    ),
):
    """Normalize standard SVCF contig names without changing coordinates.

    This command only converts standard human contig names between chr-prefixed
    and non-chr styles, e.g. 1 <-> chr1 and MT <-> chrM. Non-standard contigs
    are preserved unchanged. This is not genome-build liftover.
    """
    style = style.strip().lower()
    if style not in {"chr", "nochr"}:
        typer.echo("Error: --style must be either 'chr' or 'nochr'.", err=True)
        raise typer.Exit(code=1)

    if not dry_run and output_file is None:
        typer.echo(
            "Error: -o/--output-file is required unless --dry-run is used.",
            err=True,
        )
        raise typer.Exit(code=1)

    try:
        summary = SVCFContigNormalizer(
            input_file=input_file,
            output_file=output_file,
            style=style,
            dry_run=dry_run,
        ).run()
    except (FileNotFoundError, ValueError) as exc:
        typer.echo(f"Error: {exc}", err=True)
        raise typer.Exit(code=1) from exc

    if json_summary:
        typer.echo(summary.to_json())
    else:
        typer.echo(summary.to_text())