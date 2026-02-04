from pathlib import Path

import typer

from octopusv.formatter.svcf_to_bed_converter import SVCFtoBEDConverter
from octopusv.utils.svcf_parser import SVCFFileEventCreator


def svcf2bed(
    input_file: Path = typer.Option(..., "--input-file", "-i", exists=True, help="Input SVCF file to convert."),
    output_file: Path = typer.Option(..., "--output-file", "-o", help="Output BED file."),
    minimal: bool = typer.Option(False, "--minimal", help="Output minimal BED format (only chrom, start, end)"),
):
    """Convert SVCF file to BED format.

    The output BED file will contain structural variants with their positions.
    If --minimal is specified, only chromosome, start, and end positions will be included.
    """
    try:
        # Parse SVCF file
        typer.echo(f"Reading SVCF file: {input_file}")
        sv_event_creator = SVCFFileEventCreator([str(input_file.resolve())])
        sv_event_creator.parse()

        if not sv_event_creator.events:
            typer.echo("Warning: No events found in SVCF file", err=True)
            raise typer.Exit(code=1)

        # Convert to BED
        converter = SVCFtoBEDConverter(sv_event_creator.events, minimal=minimal)
        bed_content = converter.convert()

        # Write to output file
        with open(output_file, "w") as f:
            f.write(bed_content)

        typer.echo(f"Successfully converted SVCF to {'minimal ' if minimal else ''}BED format")
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
