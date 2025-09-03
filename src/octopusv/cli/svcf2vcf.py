from pathlib import Path

import typer

from octopusv.formatter.svcf_to_vcf_converter import SVCFtoVCFConverter
from octopusv.utils.svcf_parser import SVCFFileEventCreator


def svcf2vcf(
    input_file: Path = typer.Option(..., "--input-file", "-i", exists=True, help="Input SVCF file to convert."),
    output_file: Path = typer.Option(..., "--output-file", "-o", help="Output VCF file."),
):
    """Convert SVCF file to VCF format."""
    try:
        sv_event_creator = SVCFFileEventCreator([str(input_file.resolve())])
        sv_event_creator.parse()

        # Convert to VCF
        converter = SVCFtoVCFConverter(sv_event_creator.events, str(input_file.resolve()))
        vcf_content = converter.convert()

        # Write to output file
        with open(output_file, "w") as f:
            f.write(vcf_content)

        typer.echo(f"Converted SVCF to VCF. Output written to {output_file}")
        typer.echo("")
        typer.echo("SUCCESS: Standard VCF4.2 format generated.")
        typer.echo("This file is now compatible with bcftools, vcftools, and other standard tools.")
        typer.echo("You can proceed with downstream analysis using standard VCF workflows.")

    except FileNotFoundError:
        typer.echo(f"Error: Input file '{input_file}' not found.", err=True)
        raise typer.Exit(code=1)
    except Exception as e:
        typer.echo(f"Error occurred: {e!s}", err=True)
        raise typer.Exit(code=1)


if __name__ == "__main__":
    typer.run(svcf2vcf)
