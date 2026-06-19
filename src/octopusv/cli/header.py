"""`octopusv header` command.

Prints the header/meta contract of an SVCF/VCF file. By default it dumps every
'#' line verbatim (stopping at the first data line). With --json it emits a
structured *declared contract* parsed from the header only -- it never scans SV
records, so it does not report record counts, SV-type stats, or pass/fail
(those belong to `stat` and `validate-svcf` respectively).
"""

from pathlib import Path

import typer

from octopusv.utils.header_reader import HeaderReader


def header(
        input_svcf: Path | None = typer.Argument(
            None, exists=False, dir_okay=False, resolve_path=True,
            help="Input SVCF/VCF file."
        ),
        input_option: Path | None = typer.Option(
            None, "--input-file", "-i", dir_okay=False, resolve_path=True,
            help="Input SVCF/VCF file."
        ),
        as_json: bool = typer.Option(
            False, "--json",
            help="Emit the declared header contract as JSON instead of a raw dump."
        ),
):
    """Show the header/meta contract of an SVCF/VCF file."""
    if input_svcf and input_option:
        typer.echo(
            "Error: Specify input either as an argument or with -i/--input-file, not both.",
            err=True,
        )
        raise typer.Exit(code=1)
    input_file = input_svcf or input_option
    if input_file is None:
        typer.echo("Error: Input file is required.", err=True)
        raise typer.Exit(code=1)

    reader = HeaderReader(str(input_file))
    reader.read()

    if reader.unreadable:
        typer.echo(f"Error: Cannot read file: {input_file}", err=True)
        raise typer.Exit(code=2)

    if as_json:
        typer.echo(reader.to_json())
    else:
        typer.echo(reader.to_raw())
