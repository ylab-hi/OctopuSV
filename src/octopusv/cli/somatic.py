from pathlib import Path
import typer


def somatic(
        tumor_file: Path = typer.Option(..., "--tumor", "-t", exists=True, help="Tumor SVCF file."),
        normal_file: Path = typer.Option(..., "--normal", "-n", exists=True, help="Normal SVCF file."),
        output_file: Path = typer.Option(..., "--output-file", "-o", help="Output somatic SV file."),
        # Reuse merge parameters for SV matching
        max_distance: int | None = typer.Option(
            None,
            "--max-distance",
            help=(
                "Override the maximum breakpoint distance used to match tumor and normal SVs. "
                "If not provided, OctopuSV uses SV-type- and size-aware adaptive thresholds."
            ),
        ),
        max_length_ratio: float | None = typer.Option(
            None,
            "--max-length-ratio",
            help=(
                "Override the maximum SV length ratio used to match tumor and normal SVs. "
                "If not provided, OctopuSV uses SV-type-specific thresholds."
            ),
        ),
        min_jaccard: float = typer.Option(
            0.10,
            "--min-jaccard",
            help=(
                "Minimum interval Jaccard overlap required for DEL, DUP, and INV matching. "
                "Use 0 to disable the Jaccard requirement."
            ),
        ),
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
