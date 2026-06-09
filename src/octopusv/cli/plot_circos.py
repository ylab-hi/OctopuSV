import re
from pathlib import Path

import typer

from octopusv.utils.svcf_parser import SVCFFileEventCreator
from octopusv.vis.circos_plotter import (
    CircosPlotter,
    load_chrom_sizes_from_fai,
)


def plot_circos(
    input_file: Path = typer.Option(
        ..., "--input-file", "-i", exists=True,
        help="Input SVCF file to visualize.",
    ),
    output_file: Path = typer.Option(
        ..., "--output-file", "-o",
        help="Output figure path (.pdf / .png / .svg).",
    ),
    sample: str = typer.Option(
        None, "--sample",
        help="Sample name shown in the plot center. Defaults to input file basename.",
    ),
    fai: Path = typer.Option(
        None, "--fai",
        help="Reference FASTA .fai index for chromosome sizes. "
             "If not provided, built-in hg38 sizes are used.",
    ),

    # ---- display filters ----
    require_pass: bool = typer.Option(
        False, "--require-pass",
        help="Keep only FILTER=PASS records. Off by default, because OctopuSV "
             "merges callers whose FILTER column is not standard PASS.",
    ),
    tra_support: int = typer.Option(
        3, "--tra-support",
        help="Minimum SUPPORT for TRA (interchromosomal) links.",
    ),
    intra_support: int = typer.Option(
        3, "--intra-support",
        help="Minimum SUPPORT for intra DEL/DUP/INV links.",
    ),
    intra_min_span: int = typer.Option(
        5000, "--intra-min-span",
        help="Minimum span (bp) for intra DEL/DUP/INV links.",
    ),
    intra_max_span: int = typer.Option(
        50_000_000, "--intra-max-span",
        help="Maximum span (bp) for intra links. Larger events go to an "
             "oversized table instead of being drawn.",
    ),
    include_ins: bool = typer.Option(
        False, "--include-ins",
        help="Include INS in the DENSITY layer only (never drawn as links).",
    ),
    ins_support: int = typer.Option(
        3, "--ins-support", help="Minimum SUPPORT for INS when --include-ins is set.",
    ),
    bin_size: int = typer.Option(
        5_000_000, "--bin", help="Breakpoint density histogram bin size (bp).",
    ),

    # ---- type switches ----
    tra_only: bool = typer.Option(False, "--tra-only", help="Plot only TRA links."),
    intra_only: bool = typer.Option(False, "--intra-only", help="Plot only DEL/DUP/INV links."),
    drop_del: bool = typer.Option(False, "--drop-del", help="Drop DEL links."),
    drop_dup: bool = typer.Option(False, "--drop-dup", help="Drop DUP links."),
    drop_inv: bool = typer.Option(False, "--drop-inv", help="Drop INV links."),

    # ---- rendering ----
    tra_alpha: float = typer.Option(0.30, "--tra-alpha", help="TRA link transparency."),
    intra_alpha: float = typer.Option(0.65, "--intra-alpha", help="Intra link transparency."),
    link_lw: float = typer.Option(0.6, "--link-lw", help="Intra link line width."),
    tra_lw: float = typer.Option(1.1, "--tra-lw", help="TRA link line width (thicker)."),
    intra_height: float = typer.Option(
        0.30, "--intra-height", help="Intra arc height ratio (lower = flatter, near rim).",
    ),
    tra_height: float = typer.Option(
        0.55, "--tra-height", help="TRA arc height ratio (higher = bulges through center).",
    ),
    dpi: int = typer.Option(300, "--dpi", help="Figure DPI for raster output."),
):
    """Draw a genome-wide SV Circos overview from an SVCF file.

    Two layers are built from the SAME display-filtered SV set: an inner LINK
    layer (one arc per SV, colored by type: DEL/DUP/INV/TRA) and an outer
    breakpoint DENSITY histogram. INS is excluded by default.

    An oversized-intra table (events larger than --intra-max-span) is always
    written next to the figure for manual inspection.
    """
    # Mutually-exclusive convenience switches resolve to per-type drops.
    if tra_only and intra_only:
        typer.echo("Error: --tra-only and --intra-only are mutually exclusive.", err=True)
        raise typer.Exit(code=1)
    if tra_only:
        drop_del = drop_dup = drop_inv = True
        no_tra = False
    elif intra_only:
        no_tra = True
    else:
        no_tra = False

    # Chromosome sizes
    chrom_sizes = None
    if fai:
        chrom_sizes = load_chrom_sizes_from_fai(str(fai))
        if not chrom_sizes:
            typer.echo(f"Error: no chr1-22/X/Y sizes found in {fai}.", err=True)
            raise typer.Exit(code=1)

    sample_name = sample or re.sub(r"\.(s)?vcf(\.gz)?$", "", input_file.name)

    # Parse SVCF using OctopuSV's own parser
    creator = SVCFFileEventCreator([str(input_file)])
    creator.parse()
    if not creator.events:
        typer.echo(f"Error: no SV records parsed from {input_file}.", err=True)
        raise typer.Exit(code=1)

    plotter = CircosPlotter(
        creator.events,
        chrom_sizes=chrom_sizes,
        require_pass=require_pass,
        tra_support=tra_support,
        intra_support=intra_support,
        intra_min_span=intra_min_span,
        intra_max_span=intra_max_span,
        include_ins=include_ins,
        ins_support=ins_support,
        bin_size=bin_size,
        drop_del=drop_del,
        drop_dup=drop_dup,
        drop_inv=drop_inv,
        no_tra=no_tra,
    )
    plotter.extract()
    plotter.plot(
        str(output_file),
        sample_name=sample_name,
        tra_alpha=tra_alpha,
        intra_alpha=intra_alpha,
        link_lw=link_lw,
        tra_lw=tra_lw,
        intra_height=intra_height,
        tra_height=tra_height,
        dpi=dpi,
    )

    # Always write the oversized-intra table alongside the figure.
    out_str = str(output_file)
    oversized_path = re.sub(r"\.(pdf|png|svg)$", "", out_str) + ".oversized_intra.tsv"
    plotter.write_oversized_table(oversized_path)

    n = len(plotter.links)
    n_tra = sum(1 for r in plotter.links if r["svtype"] == "TRA")
    typer.echo(f"Circos figure written to {output_file}")
    typer.echo(f"Plotted {n} links (intra {n - n_tra}, TRA {n_tra}); "
               f"{len(plotter.breakpoints)} breakpoints in density layer.")
    typer.echo(f"Oversized intra table ({len(plotter.oversized)} events > "
               f"{intra_max_span / 1e6:g} Mb) written to {oversized_path}")


if __name__ == "__main__":
    typer.run(plot_circos)