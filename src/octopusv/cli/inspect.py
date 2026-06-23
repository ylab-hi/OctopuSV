# src/octopusv/cli/inspect.py

from __future__ import annotations

import json
from pathlib import Path
from typing import Optional

import typer

from octopusv.inspection.svcf_inspector import SVCFInspector


def _read_ids_from_file(id_file: Path) -> list[str]:
    """Read one SV ID per line from a text file.

    Blank lines and lines starting with '#' are ignored. The whole non-comment
    line is treated as the ID after stripping leading/trailing whitespace.
    """
    ids: list[str] = []

    with open(id_file, encoding="utf-8-sig") as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            ids.append(line)

    return ids


def _deduplicate_preserving_order(ids: list[str]) -> list[str]:
    """Remove duplicate IDs while preserving the first occurrence order."""
    seen: set[str] = set()
    unique_ids: list[str] = []

    for sv_id in ids:
        if sv_id in seen:
            continue
        seen.add(sv_id)
        unique_ids.append(sv_id)

    return unique_ids


def _to_jsonl_lines(result: dict) -> list[str]:
    """Render inspect results as homogeneous JSONL.

    Each output line is one inspected query result. File-level metadata is copied
    onto every line so downstream agents/tools can process JSONL without needing
    a special header/meta row with a different schema.
    """
    common = {
        "tool": result.get("tool", "octopusv inspect"),
        "input": result.get("input"),
        "file_layout": result.get("file_layout"),
        "requested_total": result.get("requested"),
        "found_total": result.get("found"),
    }

    lines: list[str] = []
    for row in result.get("results", []):
        payload = dict(common)
        payload.update(row)
        lines.append(json.dumps(payload, ensure_ascii=False, separators=(",", ":")))

    return lines


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
    id_file: Optional[Path] = typer.Option(
        None,
        "--id-file",
        exists=True,
        dir_okay=False,
        resolve_path=True,
        help=(
            "Text file containing one SV ID per line. "
            "Blank lines and lines starting with '#' are ignored."
        ),
    ),
    as_json: bool = typer.Option(
        False,
        "--json",
        help="Emit the parsed record structure as one JSON object.",
    ),
    as_jsonl: bool = typer.Option(
        False,
        "--jsonl",
        help="Emit one JSON object per inspected ID, suitable for agents and pipelines.",
    ),
):
    """Inspect one or more SVCF records by ID and render their parsed structure.

    Exposes record-level SVCF fields (endpoints, span, provenance, per-caller or
    per-sample evidence) for debugging and downstream/agent use. It does NOT
    perform gene annotation, consequence prediction, or biological interpretation.
    """
    if as_json and as_jsonl:
        typer.echo("Error: choose only one of --json or --jsonl.", err=True)
        raise typer.Exit(code=1)

    requested_ids: list[str] = []

    if sv_id:
        requested_ids.extend(sv_id)

    if id_file is not None:
        try:
            requested_ids.extend(_read_ids_from_file(id_file))
        except OSError as exc:
            typer.echo(f"Error: failed to read --id-file {id_file}: {exc}", err=True)
            raise typer.Exit(code=1) from exc

    requested_ids = _deduplicate_preserving_order(requested_ids)

    if not requested_ids:
        typer.echo("Error: at least one --id or --id-file is required.", err=True)
        raise typer.Exit(code=1)

    try:
        result = SVCFInspector(str(input_file), requested_ids).run()
    except (FileNotFoundError, ValueError) as exc:
        typer.echo(f"Error: {exc}", err=True)
        raise typer.Exit(code=1) from exc

    if as_jsonl:
        for line in _to_jsonl_lines(result):
            typer.echo(line)
    elif as_json:
        typer.echo(SVCFInspector.to_json(result))
    else:
        typer.echo(SVCFInspector.to_text(result))


if __name__ == "__main__":
    typer.run(inspect_svcf)