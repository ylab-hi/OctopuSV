from __future__ import annotations

from pathlib import Path
from typing import Optional

import typer

from octopusv.ploter.chromosome_plotter import ChromosomePlotter
from octopusv.ploter.size_plotter import SizePlotter
from octopusv.ploter.type_plotter import TypePlotter


def plot(
    input_file: Path = typer.Option(
        ...,
        "--input-file",
        "-i",
        exists=True,
        dir_okay=False,
        resolve_path=True,
        help="Input statistics file. Prefer stat.json, but legacy stat.txt is also supported.",
    ),
    output_prefix: Optional[Path] = typer.Option(
        None,
        "--output-prefix",
        "-o",
        help=(
            "Output prefix for plot files. If omitted, uses the input filename "
            "without extension."
        ),
    ),
    no_svg: bool = typer.Option(
        False,
        "--no-svg",
        help="Only write PNG files; do not write SVG files.",
    ),
):
    """Generate standard plots from OctopuSV statistics.

    Preferred input:
        octopusv stat -i sample.svcf --json -o stat.json
        octopusv plot -i stat.json -o stat

    Legacy text stat.txt is still accepted for compatibility.
    """
    prefix_path = output_prefix if output_prefix is not None else input_file.with_suffix("")
    prefix = str(prefix_path)
    save_svg = not no_svg

    ChromosomePlotter(input_file).plot(
        f"{prefix}_chromosome_distribution",
        save_svg=save_svg,
    )
    TypePlotter(input_file).plot(
        f"{prefix}_sv_types",
        save_svg=save_svg,
    )
    SizePlotter(input_file).plot(
        f"{prefix}_sv_sizes",
        save_svg=save_svg,
    )

    typer.echo(f"Plots written with prefix: {prefix}")


if __name__ == "__main__":
    typer.run(plot)
