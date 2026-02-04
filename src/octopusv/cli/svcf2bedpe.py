from pathlib import Path

import typer

from octopusv.formatter.svcf_to_bedpe_converter import SVCFtoBEDPEConverter
from octopusv.utils.svcf_parser import SVCFFileEventCreator


def svcf2bedpe(
    input_file: Path = typer.Option(..., "--input-file", "-i", exists=True, help="Input SVCF file to convert."),
    output_file: Path = typer.Option(..., "--output-file", "-o", help="Output BEDPE file."),
    minimal: bool = typer.Option(False, "--minimal", help="Output minimal BEDPE format (only coordinate columns)"),
):
    """Convert SVCF file to BEDPE format.

    The output BEDPE file will contain paired-end structural variants.
    If --minimal is specified, only coordinate columns will be included.
    """
    try:
        # Parse SVCF file
        typer.echo(f"Reading SVCF file: {input_file}")
        sv_event_creator = SVCFFileEventCreator([str(input_file.resolve())])
        sv_event_creator.parse()

        if not sv_event_creator.events:
            typer.echo("Warning: No events found in SVCF file", err=True)
            raise typer.Exit(code=1)

        # Convert to BEDPE
        converter = SVCFtoBEDPEConverter(sv_event_creator.events, minimal=minimal)
        bedpe_content = converter.convert()

        # Write to output file
        with open(output_file, "w") as f:
            f.write(bedpe_content)

        typer.echo(f"Successfully converted SVCF to {'minimal ' if minimal else ''}BEDPE format")
        typer.echo(f"Output written to {output_file}")
        typer.echo(f"Total events converted: {len(sv_event_creator.events)}")

    except FileNotFoundError:
        typer.echo(f"Error: Input file '{input_file}' not found.", err=True)
        raise typer.Exit(code=1)
    except PermissionError:
        typer.echo(f"Error: Permission denied when writing to '{output_file}'", err=True)
        raise typer.Exit(code=1)
    except Exception as e:
        typer.echo(f"Error occurred: {e!s}", err=True)
        raise typer.Exit(code=1)
