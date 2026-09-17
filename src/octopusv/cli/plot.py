from __future__ import annotations

from pathlib import Path
from typing import Optional

import typer

# Plot backends are intentionally loaded lazily so non-plot CLI commands do
# not import matplotlib or trigger font-cache initialization.  The module-level
# names remain monkeypatchable for existing tests/programmatic callers.
ChromosomePlotter = None
TypePlotter = None
SizePlotter = None


def _load_plotters():
    global ChromosomePlotter, TypePlotter, SizePlotter
    if ChromosomePlotter is None or TypePlotter is None or SizePlotter is None:
        from octopusv.ploter.chromosome_plotter import ChromosomePlotter as _ChromosomePlotter
        from octopusv.ploter.size_plotter import SizePlotter as _SizePlotter
        from octopusv.ploter.type_plotter import TypePlotter as _TypePlotter

        ChromosomePlotter = _ChromosomePlotter
        TypePlotter = _TypePlotter
        SizePlotter = _SizePlotter
    return ChromosomePlotter, TypePlotter, SizePlotter


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
    chromosome_plotter, type_plotter, size_plotter = _load_plotters()

    chromosome_plotter(input_file).plot(
        f"{prefix}_chromosome_distribution",
        save_svg=save_svg,
    )
    type_plotter(input_file).plot(
        f"{prefix}_sv_types",
        save_svg=save_svg,
    )
    size_plotter(input_file).plot(
        f"{prefix}_sv_sizes",
        save_svg=save_svg,
    )

    typer.echo(f"Plots written with prefix: {prefix}")


if __name__ == "__main__":
    typer.run(plot)
