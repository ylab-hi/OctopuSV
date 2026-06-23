# src/octopusv/cli/inspect.py

from __future__ import annotations

from pathlib import Path
from typing import Optional

import typer

from octopusv.inspection.svcf_inspector import SVCFInspector


def inspect_svcf(
    input_file: Path = typer.Option(
        ...,
        "--input-file",
        "-i",
        exists=True,
        dir_okay=False,
        resolve_path=True,
        help="Input SVCF file.",
    ),
    sv_id: Optional[list[str]] = typer.Option(
        None,
        "--id",
        help="SV ID to inspect. Can be used multiple times.",
    ),
    as_json: bool = typer.Option(
        False,
        "--json",
        help="Emit the parsed record structure as JSON (the agent-facing path).",
    ),
):
    """Inspect one or more SVCF records by ID and render their parsed structure.

    Exposes record-level SVCF fields (endpoints, span, provenance, per-caller or
    per-sample evidence) for debugging and downstream/agent use. It does NOT
    perform gene annotation, consequence prediction, or biological interpretation.
    """
    if not sv_id:
        typer.echo("Error: at least one --id is required.", err=True)
        raise typer.Exit(code=1)

    try:
        result = SVCFInspector(str(input_file), sv_id).run()
    except (FileNotFoundError, ValueError) as exc:
        typer.echo(f"Error: {exc}", err=True)
        raise typer.Exit(code=1) from exc

    if as_json:
        typer.echo(SVCFInspector.to_json(result))
    else:
        typer.echo(SVCFInspector.to_text(result))


if __name__ == "__main__":
    typer.run(inspect_svcf)