import base64
import logging
import sys
from datetime import datetime
from pathlib import Path

from jinja2 import Template
from octopusv import __PACKAGE_NAME__

logging.basicConfig(level=logging.INFO)

def load_template() -> Path:
    """Load HTML template for report generation."""
    template_path = Path(sys.modules[__PACKAGE_NAME__].__file__)
    if template_path.parent is not None:
        template_path = template_path.parent / "report"
        return template_path / "template.html"
    raise FileNotFoundError("Could not locate template directory")

def load_logo() -> Path:
    """Load logo for report generation."""
    logo_path = Path(sys.modules[__PACKAGE_NAME__].__file__)
    if logo_path.parent is not None:
        logo_path = logo_path.parent / "report"
        return logo_path / "logo.png"
    raise FileNotFoundError("Could not locate logo directory")

def image_to_base64(image_path: Path) -> str:
    """Convert image to base64 string with error handling."""
    try:
        if not image_path.exists():
            logging.warning(f"Image file not found: {image_path}")
            return ""
        with image_path.open("rb") as f:
            return base64.b64encode(f.read()).decode("utf-8")
    except Exception as e:
        logging.error(f"Error converting {image_path} to base64: {str(e)}")
        return ""

class ReportGenerator:
    """Generate HTML report using Jinja2 template."""
    def __init__(self, template_path=None):
        self.template_path = template_path or load_template()
        with Path(self.template_path).open() as f:
            self.template = Template(f.read())

    def _get_plot_images(self, base_path: str) -> dict:
        """Get base64 encoded images from existing plot files."""
        plots = {}
        plot_files = {
            'chromosome_coverage_plot_base64': f"{base_path}_chromosome_distribution.png",
            'sv_distribution_plot_base64': f"{base_path}_sv_types.png",
            'size_distribution_plot_base64': f"{base_path}_sv_sizes.png"
        }
        
        for key, filepath in plot_files.items():
            try:
                plots[key] = image_to_base64(Path(filepath))
            except Exception as e:
                logging.error(f"Error processing {key}: {str(e)}")
        
        return plots

    def generate(self, input_file: str, output_path: str, sample_id: str, summary_stats: dict):
        """Generate HTML report with integrated plots and statistics."""
        try:
            # Get plot images
            output_prefix = str(Path(output_path).with_suffix(''))
            plots = self._get_plot_images(output_prefix)
            
            # Prepare template data
            template_data = {
                "generation_date": datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
                "sample_id": sample_id,
                "logo_path": load_logo(),
                **plots,
                **summary_stats
            }

            # Render template and write output
            html_content = self.template.render(**template_data)
            with Path(output_path).with_suffix('.html').open("w", encoding='utf-8') as f:
                f.write(html_content)
            logging.info(f"Successfully generated HTML report at {output_path}")

        except Exception as e:
            logging.error(f"Error generating HTML report: {str(e)}")
            raise
