# src/octopusv/cli/subset.py

from __future__ import annotations

from pathlib import Path
from typing import Optional

import typer

from octopusv.subset.svcf_subset import (
    SVCFSubset,
    SubsetConfig,
    expand_csv_values,
    read_name_file,
    unique_preserve,
)


def subset_svcf(
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
        help="Output subset SVCF file. Required unless --dry-run is used.",
    ),
    sample: Optional[list[str]] = typer.Option(
        None,
        "--sample",
        help=(
            "Sample/input column to retain in sample-mode SVCF. "
            "Can be used multiple times or as comma-separated values."
        ),
    ),
    sample_file: Optional[Path] = typer.Option(
        None,
        "--sample-file",
        exists=True,
        dir_okay=False,
        resolve_path=True,
        help="File containing sample/input column names to retain, one per line.",
    ),
    caller: Optional[list[str]] = typer.Option(
        None,
        "--caller",
        help=(
            "Caller evidence block to retain in caller-mode SVCF. "
            "Can be used multiple times or as comma-separated values."
        ),
    ),
    caller_file: Optional[Path] = typer.Option(
        None,
        "--caller-file",
        exists=True,
        dir_okay=False,
        resolve_path=True,
        help="File containing caller names to retain, one per line.",
    ),
    mode: str = typer.Option(
        "auto",
        "--mode",
        help="Subset mode: auto, sample, or caller.",
    ),
    keep_empty: bool = typer.Option(
        False,
        "--keep-empty",
        help=(
            "Keep records not supported by selected samples/callers. "
            "Default behavior is to drop empty records."
        ),
    ),
    no_audit_header: bool = typer.Option(
        False,
        "--no-audit-header",
        help="Do not add OctopuSV subset audit header lines.",
    ),
    json_summary: bool = typer.Option(
        False,
        "--json-summary",
        help="Print an agent-friendly JSON subset summary.",
    ),
    dry_run: bool = typer.Option(
        False,
        "--dry-run",
        help="Do not write output SVCF; only report subset summary.",
    ),
):
    """Subset SVCF sample/caller evidence columns while preserving SVCF structure."""
    if not dry_run and output_file is None:
        typer.echo(
            "Error: -o/--output-file is required unless --dry-run is used.",
            err=True,
        )
        raise typer.Exit(code=1)

    selected_samples = unique_preserve(
        expand_csv_values(sample) + read_name_file(sample_file)
    )
    selected_callers = unique_preserve(
        expand_csv_values(caller) + read_name_file(caller_file)
    )

    if selected_samples and selected_callers:
        typer.echo(
            "Error: Use either --sample/--sample-file or --caller/--caller-file, not both.",
            err=True,
        )
        raise typer.Exit(code=1)

    if not selected_samples and not selected_callers:
        typer.echo(
            "Error: At least one --sample or --caller is required.",
            err=True,
        )
        raise typer.Exit(code=1)

    try:
        config = SubsetConfig(
            input_file=str(input_file),
            output_file=str(output_file) if output_file is not None else None,
            dry_run=dry_run,
            mode=mode,
            selected_samples=selected_samples,
            selected_callers=selected_callers,
            keep_empty=keep_empty,
            audit_header=not no_audit_header,
        )

        summary = SVCFSubset(config).run()

    except (FileNotFoundError, ValueError) as exc:
        typer.echo(f"Error: {exc}", err=True)
        raise typer.Exit(code=1) from exc

    if json_summary:
        typer.echo(SVCFSubset.summary_to_json(summary))
    else:
        typer.echo(SVCFSubset.summary_to_text(summary))


if __name__ == "__main__":
    typer.run(subset_svcf)