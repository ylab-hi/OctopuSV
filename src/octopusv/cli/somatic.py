from pathlib import Path
import typer


def somatic(
        tumor_file: Path = typer.Option(..., "--tumor", "-t", exists=True, help="Tumor SVCF file."),
        normal_file: Path = typer.Option(..., "--normal", "-n", exists=True, help="Normal SVCF file."),
        output_file: Path = typer.Option(..., "--output-file", "-o", help="Output somatic SV file."),

        # Reuse merge parameters for SV matching
        max_distance: int = typer.Option(50, "--max-distance", help="Maximum allowed distance for merging events."),
        max_length_ratio: float = typer.Option(1.3, "--max-length-ratio",
                                               help="Maximum allowed ratio between event lengths."),
        min_jaccard: float = typer.Option(0.7, "--min-jaccard", help="Minimum required Jaccard index for overlap."),
):
    """Extract somatic SVs by finding tumor-specific variants (tumor - normal intersection)."""

    # Import merge function and call it directly
    from octopusv.cli.merge import merge

    # Call merge with specific logic: extract SVs only present in tumor
    merge(
        input_files=[tumor_file, normal_file],
        input_option=None,  # Set to None explicitly
        output_file=output_file,
        mode="sample",
        sample_names="tumor,normal",
        specific=[tumor_file],  # Extract tumor-specific variants only
        max_distance=max_distance,
        max_length_ratio=max_length_ratio,
        min_jaccard=min_jaccard,
        intersect=False,
        union=False,
        upsetr=False,
        min_support=None,
        exact_support=None,
        max_support=None,
        expression=None,
        tra_delta=50,
        tra_min_overlap_ratio=0.5,
        tra_strand_consistency=True,
        bnd_delta=50,
        caller_names=None,
        upsetr_output=None
    )

    typer.echo(f"Somatic SVs extracted and written to: {output_file}")
    typer.echo("Note: Output contains SVs found in tumor but not in normal sample.")
    typer.echo("Remember to convert SVCF to VCF format for downstream analysis:")
    typer.echo(f"  octopusv svcf2vcf -i {output_file} -o {output_file.with_suffix('.vcf')}")