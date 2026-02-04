from pathlib import Path
import typer
from octopusv.report.generator import ReportGenerator
from octopusv.stater.sv_stater import SVStater
from octopusv.ploter.chromosome_plotter import ChromosomePlotter
from octopusv.ploter.type_plotter import TypePlotter
from octopusv.ploter.size_plotter import SizePlotter

def stat(
    input_file: Path = typer.Argument(..., help="Input SVCF file to analyze."),
    output_file: Path = typer.Option(..., "--output-file", "-o", help="Output file for statistics."),
    min_size: int = typer.Option(50, "--min-size", help="Minimum SV size to consider."),
    max_size: int = typer.Option(None, "--max-size", help="Maximum SV size to consider."),
    report: bool = typer.Option(False, "--report", help="Generate an HTML report."),
):
    """Analyze a single SVCF file and generate comprehensive statistics."""
    # Run analysis
    sv_stater = SVStater(str(input_file), min_size=min_size, max_size=max_size)
    sv_stater.analyze()
    sv_stater.write_results(output_file)

    if report:
        typer.echo("Generating HTML report...")
        
        # Generate plots
        output_prefix = str(output_file.with_suffix(''))
        
        # Generate chromosome distribution plot
        chromosome_plotter = ChromosomePlotter(str(output_file))
        chromosome_plotter.plot(f"{output_prefix}_chromosome_distribution")
        
        # Generate SV type plot
        type_plotter = TypePlotter(str(output_file))
        type_plotter.plot(f"{output_prefix}_sv_types")
        
        # Generate size distribution plot
        size_plotter = SizePlotter(str(output_file))
        size_plotter.plot(f"{output_prefix}_sv_sizes")
        
        # Generate HTML report
        report_generator = ReportGenerator()
        summary_stats = sv_stater.export_html(output_file)
        report_generator.generate(
            input_file=str(input_file),
            output_path=str(output_file),
            sample_id=input_file.stem,
            summary_stats=summary_stats
        )
        typer.echo("Report generated.")

    typer.echo(f"Analysis results written to {output_file}")

if __name__ == "__main__":
    typer.run(stat)
