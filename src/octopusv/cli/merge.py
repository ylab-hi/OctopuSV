from pathlib import Path

import typer

from octopusv.merger.sv_merger import SVMerger
from octopusv.merger.upset_plotter import UpSetPlotter  # Import UpSetPlotter
from octopusv.utils.SV_classifier_by_chromosome import SVClassifiedByChromosome
from octopusv.utils.SV_classifier_by_type import SVClassifierByType
from octopusv.utils.svcf_parser import SVCFFileEventCreator


def get_contigs_from_svcf(filenames):
    """Extract contig information from SVCF files.

    Args:
        filenames (list): List of SVCF filenames.

    Returns:
        dict: Dictionary of contig IDs and their lengths.
    """
    contigs = {}
    for filename in filenames:
        with open(filename) as f:
            for line in f:
                if line.startswith("##contig"):
                    # Parse contig information
                    line = line.strip()
                    if line.startswith("##contig=<") and line.endswith(">"):
                        content = line[len("##contig=<") : -1]
                        parts = content.split(",")
                        contig_id = ""
                        contig_length = ""
                        for part in parts:
                            if part.startswith("ID="):
                                contig_id = part.split("=", 1)[1]
                            elif part.startswith("length="):
                                contig_length = part.split("=", 1)[1]
                        if contig_id and contig_length:
                            contigs[contig_id] = contig_length
                elif not line.startswith("##"):
                    break
    return contigs


def merge(
    input_files: list[Path] = typer.Argument(None, help="List of input SVCF files to merge."),
    input_option: list[Path] = typer.Option(None, "--input-file", "-i", help="Input SVCF files to merge."),
    output_file: Path = typer.Option(..., "--output-file", "-o", help="Output file for merged SV data."),
    intersect: bool = typer.Option(False, "--intersect", help="Apply intersection strategy for merging."),
    union: bool = typer.Option(False, "--union", help="Apply union strategy for merging."),
    specific: list[Path] = typer.Option(
        None, "--specific", help="Extract SVs that are specifically supported by provided files."
    ),
    min_support: int = typer.Option(None, "--min-support", help="Minimum number of files that must support an SV."),
    exact_support: int = typer.Option(None, "--exact-support", help="Exact number of files that must support an SV."),
    max_support: int = typer.Option(None, "--max-support", help="Maximum number of files that can support an SV."),
    expression: str = typer.Option(
        None,
        "--expression",
        help="Logical expression for complex file combinations (e.g., '(A AND B) AND NOT (C OR D)')",
    ),
    max_distance: int = typer.Option(
        50, "--max-distance", help="Maximum allowed distance between start or end positions for merging events."
    ),
    max_length_ratio: float = typer.Option(
        1.3, "--max-length-ratio", help="Maximum allowed ratio between event lengths for merging events."
    ),
    min_jaccard: float = typer.Option(
        0.7, "--min-jaccard", help="Minimum required Jaccard index for overlap to merge events."
    ),
    tra_delta: int = typer.Option(
        50, "--tra-delta", help="Position uncertainty threshold for TRA events (in base pairs)."
    ),
    tra_min_overlap_ratio: float = typer.Option(0.5, "--tra-min-overlap", help="Minimum overlap ratio for TRA events."),
    tra_strand_consistency: bool = typer.Option(
        True, "--tra-strand-consistency", help="Whether to require strand consistency for TRA events."
    ),
    bnd_delta: int = typer.Option(
        50, "--bnd-delta", help="Position uncertainty threshold for BND events (in base pairs)."
    ),
    upsetr: bool = typer.Option(
        False, "--upsetr", help="Generate UpSet plot visualization of input file intersections."
    ),
    upsetr_output: Path = typer.Option(
        None,
        "--upsetr-output",
        help="Output path for UpSet plot. If not provided, will use output_file basename with _upset.png suffix.",
    ),

):
    """Merge multiple SVCF files based on specified strategy."""
    # Combine input files from arguments and options
    all_input_files = (input_files or []) + (input_option or [])

    if not all_input_files:
        typer.echo("Error: No input files provided.", err=True)
        raise typer.Exit(code=1)
    if specific and not specific[0]:
        typer.echo("Error: --specific option requires at least one file.", err=True)
        raise typer.Exit(code=1)
    if min_support is not None and min_support < 1:
        typer.echo("Error: --min-support must be a positive integer.", err=True)
        raise typer.Exit(code=1)

    # Get contig information from input files
    input_filenames = [str(file) for file in all_input_files]
    contigs = get_contigs_from_svcf(input_filenames)

    sv_event_creator = SVCFFileEventCreator(input_filenames)
    sv_event_creator.parse()
    classifier = SVClassifierByType(sv_event_creator.events)
    classifier.classify()
    chromosome_classifier = SVClassifiedByChromosome(classifier.get_classified_events())
    chromosome_classifier.classify()

    sv_merger = SVMerger(
        chromosome_classifier.get_classified_events(),
        all_input_files=input_filenames,
        tra_delta=tra_delta,
        tra_min_overlap_ratio=tra_min_overlap_ratio,
        tra_strand_consistency=tra_strand_consistency,
        max_distance=max_distance,
        max_length_ratio=max_length_ratio,
        min_jaccard=min_jaccard,
        bnd_delta=bnd_delta,
    )
    sv_merger.merge()

    if expression:
        results = sv_merger.get_events_by_expression(expression)
    elif intersect:
        results = sv_merger.get_events_by_source([str(file) for file in all_input_files], operation="intersection")
    elif union:
        results = sv_merger.get_events_by_source([str(file) for file in all_input_files], operation="union")
    elif specific:
        specific_files = [str(file) for file in specific]
        results = sv_merger.get_events_by_source(specific_files, operation="specific")
    elif exact_support is not None:
        results = sv_merger.get_events_by_exact_support(exact_support)
    elif min_support is not None or max_support is not None:
        results = sv_merger.get_events_by_support_range(min_support, max_support)
    else:
        raise ValueError(
            "No merge strategy specified. Please use --intersect, --union, --specific, "
            "--min-support, --exact-support, --max-support, or --expression."
        )

    # Write results with contig information
    sv_merger.write_results(output_file, results, contigs)
    typer.echo(f"Merged results written to {output_file}")

    # Generate UpSet plot if requested
    if upsetr:
        try:
            plot_file = str(upsetr_output) if upsetr_output else str(output_file).rsplit(".", 1)[0] + "_upset.png"
            plotter = UpSetPlotter(sv_merger.get_all_merged_events(), all_input_files)
            plotter.plot(plot_file)
            typer.echo(f"UpSet plot written to {plot_file}")
        except ImportError:
            typer.echo(
                "Warning: Could not generate UpSet plot. Please ensure matplotlib and numpy are installed.", err=True
            )
        except Exception as e:
            typer.echo(f"Warning: Failed to generate UpSet plot: {e!s}", err=True)


if __name__ == "__main__":
    typer.run(merge)
