import logging
from pathlib import Path

import typer
from rich.logging import RichHandler

from octopusv.bencher.sv_bencher import SVBencher

# Set logging style
FORMAT = "%(message)s"
logging.basicConfig(level="INFO", format=FORMAT, datefmt="[%X]", handlers=[RichHandler()])


def bench(
    truth_file: Path = typer.Argument(..., help="Path to the truth (ground truth) SVCF file."),
    call_file: Path = typer.Argument(..., help="Path to the call (test) SVCF file."),
    output_dir: Path = typer.Option(..., "--output-dir", "-o", help="Output directory for benchmark results."),
    reference_distance: int = typer.Option(
        500, "--reference-distance", "-r", help="Max reference location distance (default: 500)"
    ),
    size_similarity: float = typer.Option(
        0.7, "--size-similarity", "-P", help="Min pct size similarity (minsize/maxsize) (default: 0.7)"
    ),
    reciprocal_overlap: float = typer.Option(
        0.0, "--reciprocal-overlap", "-O", help="Min pct reciprocal overlap (default: 0.0)"
    ),
    type_ignore: bool = typer.Option(
        False, "--type-ignore", "-t", help="Variant types don't need to match to compare (default: False)"
    ),
    size_min: int = typer.Option(
        50, "--size-min", "-s", help="Minimum variant size to consider from test calls (default: 50)"
    ),
    size_max: int = typer.Option(50000, "--size-max", help="Maximum variant size to consider (default: 50000)"),
    pass_only: bool = typer.Option(
        False, "--pass-only", help="Only consider variants with FILTER == PASS (default: False)"
    ),
    enable_sequence_comparison: bool = typer.Option(
        False,
        "--enable-sequence-comparison",
        help="Enable sequence similarity comparison if sequences available (default: False)",
    ),
    sequence_similarity: float = typer.Option(
        0.7,
        "--sequence-similarity",
        "-p",
        help="Min sequence similarity when sequence comparison enabled (default: 0.7)",
    ),
):
    """Benchmark structural variation calls against a truth set using GIAB standards.

    This tool implements the GIAB (Genome in a Bottle) consortium's recommendations for
    structural variant benchmarking. It compares a test VCF file against a truth set
    and calculates precision, recall, and F1 scores.

    Default thresholds are set according to GIAB standards:
    - 500bp maximum reference distance
    - 70% size similarity
    - 0% reciprocal overlap
    - 50bp minimum variant size
    - 50kb maximum variant size

    Results are written to the specified output directory, including:
    - True positives from the baseline (tp-base.vcf)
    - True positives from the test set (tp-call.vcf)
    - False positives (fp.vcf)
    - False negatives (fn.vcf)
    - Summary metrics in JSON format (summary.json)
    """
    try:
        bencher = SVBencher(
            truth_file,
            call_file,
            output_dir,
            reference_distance=reference_distance,
            sequence_similarity=sequence_similarity,
            size_similarity=size_similarity,
            reciprocal_overlap=reciprocal_overlap,
            type_ignore=type_ignore,
            size_min=size_min,
            size_max=size_max,
            pass_only=pass_only,
            enable_sequence_comparison=enable_sequence_comparison,
        )

        typer.echo("Starting benchmark...")
        bencher.run_benchmark()
        typer.echo(f"Benchmark results written to {output_dir}")

    except Exception as e:
        typer.echo(f"Error during benchmarking: {e!s}", err=True)
        raise typer.Exit(code=1)
