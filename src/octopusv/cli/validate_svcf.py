"""`octopusv validate-svcf` command.

Checks whether an SVCF file conforms to the OctopuSV SVCF contract.

Legal modes:
  1. single-caller correct output
  2. caller-merge output
  3. sample/multi output

Default output is a human-readable summary.
Use --json for an agent-friendly structured report.

Exit codes:
  0 = pass / pass-with-warnings
  1 = validation failed
  2 = unreadable file
"""

from pathlib import Path

import typer

from octopusv.utils.svcf_validator import SVCFValidator


def validate_svcf(
    input_svcf: Path | None = typer.Argument(
        None,
        exists=False,
        dir_okay=False,
        resolve_path=True,
        help="Input SVCF file to validate.",
    ),
    input_option: Path | None = typer.Option(
        None,
        "--input-file",
        "-i",
        dir_okay=False,
        resolve_path=True,
        help="Input SVCF file to validate.",
    ),
    as_json: bool = typer.Option(
        False,
        "--json",
        help="Emit a structured JSON report instead of text.",
    ),
    strict_co: bool = typer.Option(
        False,
        "--strict-co",
        help="Treat unparseable CO as an error instead of a warning.",
    ),
    max_issues: int = typer.Option(
        50,
        "--max-issues",
        help=(
            "Maximum number of issues to print or include in JSON. "
            "Use 0 to show all issues."
        ),
    ),
):
    """Validate an SVCF file against the OctopuSV SVCF contract."""
    if input_svcf and input_option:
        typer.echo(
            "Error: Specify input either as an argument or with -i/--input-file, not both.",
            err=True,
        )
        raise typer.Exit(code=1)

    input_file = input_svcf or input_option
    if input_file is None:
        typer.echo("Error: Input SVCF file is required.", err=True)
        raise typer.Exit(code=1)

    if max_issues < 0:
        typer.echo("Error: --max-issues must be >= 0.", err=True)
        raise typer.Exit(code=1)

    issue_limit = None if max_issues == 0 else max_issues

    validator = SVCFValidator(str(input_file), strict_co=strict_co)
    validator.validate()

    if as_json:
        typer.echo(validator.to_json(max_issues=issue_limit))
    else:
        typer.echo(validator.to_summary(max_issues=issue_limit))

    raise typer.Exit(code=validator.exit_code())
