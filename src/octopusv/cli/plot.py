from pathlib import Path

import typer

from octopusv.ploter.chromosome_plotter import ChromosomePlotter
from octopusv.ploter.size_plotter import SizePlotter
from octopusv.ploter.type_plotter import TypePlotter


def plot(
    input_file: Path = typer.Argument(..., help="Input stat.txt file to plot."),
    output_prefix: Path = typer.Option(..., "--output-prefix", "-o", help="Output prefix for plot files."),
):
    """Generate plots from the statistics file."""
    chromosome_plotter = ChromosomePlotter(input_file)
    type_plotter = TypePlotter(input_file)
    size_plotter = SizePlotter(input_file)

    chromosome_plotter.plot(f"{output_prefix}_chromosome_distribution")
    type_plotter.plot(f"{output_prefix}_sv_types")
    size_plotter.plot(f"{output_prefix}_sv_sizes")

    typer.echo(f"Plots have been generated with prefix: {output_prefix}")


if __name__ == "__main__":
    typer.run(plot)
