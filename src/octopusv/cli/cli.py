import logging
import sys
import typer
from rich.logging import RichHandler
from octopusv import __version__
from .bench import bench
from .convert import correct
from .merge import merge
from .plot import plot
from .stat import stat
from .svcf2bed import svcf2bed
from .svcf2bedpe import svcf2bedpe
from .svcf2vcf import svcf2vcf
from .somatic import somatic

app = typer.Typer(
    epilog=f"{typer.style('Agent Octopus Code V helps you dive deep into the structural variations ocean!', fg=typer.colors.GREEN, bold=True)}",
    context_settings={"help_option_names": ["-h", "--help"]},
)

FORMAT = "%(message)s"
logging.basicConfig(level="INFO", format=FORMAT, datefmt="[%X]", handlers=[RichHandler()])

app.command()(correct)  # Command to initiate convert functionality.
app.command()(merge)  # Command to initiate merge functionality.
app.command(name="benchmark")(bench)  # Command to initiate bench functionality.
app.command()(stat)  # Command to initiate stat functionality.
app.command()(plot)  # Command to initiate plot functionality.
app.command()(somatic)  # Command to initiate somatic SV calling functionality.
app.command()(svcf2vcf)
app.command()(svcf2bed)
app.command()(svcf2bedpe)

@app.callback()
def display_version_info():
    """Display the version information."""

display_version_info.__doc__ = f"""{typer.style("Version", fg=typer.colors.YELLOW, bold=True)}: {typer.style(f"{__version__}", fg=typer.colors.GREEN, bold=True)}"""

# For documentation purposes
if "sphinx" in sys.modules and __name__ != "__main__":
    # Create the typer click object to generate docs with sphinx-click
    typer_click_object = typer.main.get_command(app)

if __name__ == "__main__":
    app()