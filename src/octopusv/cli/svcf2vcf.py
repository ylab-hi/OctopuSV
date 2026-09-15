from pathlib import Path

import typer

from octopusv.formatter.svcf_to_vcf_converter import SVCFtoVCFConverter


def svcf2vcf(
    input_file: Path = typer.Option(
        ...,
        "--input-file",
        "-i",
        exists=True,
        help="Input SVCF file to convert.",
    ),
    output_file: Path = typer.Option(
        ...,
        "--output-file",
        "-o",
        help="Output VCF file.",
    ),
):
    """Convert SVCF file to VCF format."""
    try:
        resolved_input = str(input_file.resolve())

        # Stream one SVCF record at a time.  convert_to_file writes through a
        # temporary file and atomically replaces the requested output only
        # after the entire conversion succeeds.
        converter = SVCFtoVCFConverter(
            events=None,
            input_svcf_file=resolved_input,
        )
        converter.convert_to_file(output_file)

        typer.echo(f"Converted SVCF to VCF. Output written to {output_file}")
        typer.echo("")
        typer.echo("SUCCESS: VCF4.2-compatible output generated.")
        typer.echo(
            "This file is now compatible with bcftools, vcftools, "
            "and other standard tools."
        )
        typer.echo(
            "You can proceed with downstream analysis using standard "
            "VCF workflows."
        )

    except FileNotFoundError:
        typer.echo(f"Error: Input file '{input_file}' not found.", err=True)
        raise typer.Exit(code=1)
    except Exception as exc:
        typer.echo(f"Error occurred: {exc!s}", err=True)
        raise typer.Exit(code=1)


if __name__ == "__main__":
    typer.run(svcf2vcf)
