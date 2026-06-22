import logging
import sys

import typer
from rich.logging import RichHandler

from octopusv import __version__
from .bench import bench
from .clean import clean
from .convert import correct
from .header import header
from .merge import merge
from .plot import plot
from .plot_circos import plot_circos
from .somatic import somatic
from .stat import stat
from .svcf2bed import svcf2bed
from .svcf2bedpe import svcf2bedpe
from .svcf2vcf import svcf2vcf
from .validate_svcf import validate_svcf
from .filter import filter_svcf
from .query import query_svcf
from .subset import subset_svcf


HELP = (
    f"[bold yellow]Version[/bold yellow]: [bold green]{__version__}[/bold green]\n\n"
    "OctopuSV: SVCF-centered structural variant standardization, inspection, "
    "statistics, visualization, and downstream-ready operations."
)

EPILOG = (
    "[bold green]OctopuSV[/bold green] — "
    "navigate the structural-variation ocean."
)


app = typer.Typer(
    help=HELP,
    epilog=EPILOG,
    rich_markup_mode="rich",
    context_settings={"help_option_names": ["-h", "--help"]},
)

FORMAT = "%(message)s"
logging.basicConfig(
    level="INFO",
    format=FORMAT,
    datefmt="[%X]",
    handlers=[RichHandler()],
)


# ---------------------------------------------------------------------
# SVCF inspection and statistics
# ---------------------------------------------------------------------

app.command(
    rich_help_panel="SVCF inspection and statistics",
)(header)

app.command(
    name="validate-svcf",
    rich_help_panel="SVCF inspection and statistics",
)(validate_svcf)

app.command(
    rich_help_panel="SVCF inspection and statistics",
)(stat)

app.command(
    name="filter",
    rich_help_panel="SVCF filtering and querying",
)(filter_svcf)

app.command(
    name="query",
    rich_help_panel="SVCF filtering and querying",
)(query_svcf)

app.command(
    name="subset",
    rich_help_panel="SVCF filtering and querying",
)(subset_svcf)

# ---------------------------------------------------------------------
# SV standardization and merging
# ---------------------------------------------------------------------

app.command(
    rich_help_panel="SV standardization and merging",
)(clean)

app.command(
    rich_help_panel="SV standardization and merging",
)(correct)

app.command(
    rich_help_panel="SV standardization and merging",
)(merge)


# ---------------------------------------------------------------------
# Conversion
# ---------------------------------------------------------------------

app.command(
    rich_help_panel="Conversion",
)(svcf2vcf)

app.command(
    rich_help_panel="Conversion",
)(svcf2bed)

app.command(
    rich_help_panel="Conversion",
)(svcf2bedpe)


# ---------------------------------------------------------------------
# Analysis and visualization
# ---------------------------------------------------------------------

app.command(
    name="benchmark",
    rich_help_panel="Analysis and visualization",
)(bench)

app.command(
    rich_help_panel="Analysis and visualization",
)(plot)

app.command(
    name="plot-circos",
    rich_help_panel="Analysis and visualization",
)(plot_circos)

app.command(
    rich_help_panel="Analysis and visualization",
)(somatic)


@app.callback()
def display_version_info():
    """OctopuSV command-line interface."""


# For documentation purposes
if "sphinx" in sys.modules and __name__ != "__main__":
    # Create the typer click object to generate docs with sphinx-click.
    typer_click_object = typer.main.get_command(app)


if __name__ == "__main__":
    app()
