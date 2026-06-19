from pathlib import Path

import typer

from octopusv.report.generator import ReportGenerator
from octopusv.stater.sv_stater import SVStater
from octopusv.ploter.chromosome_plotter import ChromosomePlotter
from octopusv.ploter.type_plotter import TypePlotter
from octopusv.ploter.size_plotter import SizePlotter


def stat(
    input_arg: Path = typer.Argument(
        None, help="Input SVCF file to analyze (positional). Alternatively use -i/--input-file."
    ),
    input_option: Path = typer.Option(
        None, "--input-file", "-i", help="Input SVCF file to analyze."
    ),
    output_file: Path = typer.Option(
        None, "--output-file", "-o",
        help="Output file. Optional for text/JSON (defaults to stdout); required with --report."
    ),
    as_json: bool = typer.Option(
        False, "--json",
        help="Emit structured JSON. Writes to -o if given, else to stdout."
    ),
    fai: Path = typer.Option(
        None, "--fai",
        help="Reference .fai index; highest-priority source for chromosome density."
    ),
    genome: str = typer.Option(
        "auto", "--genome",
        help="Built-in genome for density when no --fai: auto|hg38|grch38|hg19|grch37. "
             "auto detects from contig names (hs37d5->hg19); ambiguous -> hg38 (assumed)."
    ),
    min_size: int = typer.Option(50, "--min-size", help="Minimum SV size to consider."),
    max_size: int = typer.Option(None, "--max-size", help="Maximum SV size to consider."),
    report: bool = typer.Option(False, "--report", help="Generate an HTML report (requires -o)."),
):
    """Analyze a single SVCF file and report SV-content statistics."""
    # Resolve input.
    input_file = input_option if input_option is not None else input_arg
    if input_file is None:
        typer.echo(
            "Error: no input SVCF provided. Pass it as a positional argument "
            "or with -i/--input-file.",
            err=True,
        )
        raise typer.Exit(code=1)

    # --report needs an output prefix to write plots + html next to.
    if report and output_file is None:
        typer.echo("Error: --report requires -o/--output-file.", err=True)
        raise typer.Exit(code=1)

    sv_stater = SVStater(str(input_file), min_size=min_size, max_size=max_size,
                         fai=str(fai) if fai else None, genome=genome)
    sv_stater.analyze()

    # JSON path: to file if -o given, else stdout.
    if as_json:
        payload = sv_stater.write_json(str(output_file) if output_file else None)
        if payload is not None:
            typer.echo(payload)
        else:
            typer.echo(f"JSON statistics written to {output_file}")
        return

    # Text path: to file if -o given, else stdout.
    if output_file is None:
        typer.echo(sv_stater.to_text())
        return

    sv_stater.write_results(output_file)

    if report:
        typer.echo("Generating HTML report...")
        output_prefix = str(output_file.with_suffix(''))
        ChromosomePlotter(str(output_file)).plot(f"{output_prefix}_chromosome_distribution")
        TypePlotter(str(output_file)).plot(f"{output_prefix}_sv_types")
        SizePlotter(str(output_file)).plot(f"{output_prefix}_sv_sizes")
        report_generator = ReportGenerator()
        report_generator.generate(
            input_file=str(input_file),
            output_path=str(output_file),
            sample_id=input_file.stem,
            summary_stats=sv_stater.to_html_dict(),
        )
        typer.echo("Report generated.")

    typer.echo(f"Analysis results written to {output_file}")


if __name__ == "__main__":
    typer.run(stat)
