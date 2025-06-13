import logging
from pathlib import Path

import typer

from octopusv.converter.base import get_alt_chrom_pos
from octopusv.converter.bnd_keeping import BNDKeepingConverter
from octopusv.converter.mpi2tra import MatePairIndependentToTRAConverter
from octopusv.converter.mpm2tra import MatePairMergeToTRAConverter
from octopusv.converter.mprtra2tra import MatePairReciprocalTranslocationToTRAConverter
from octopusv.converter.nobnd import NonBNDConverter
from octopusv.converter.snmd_dndpi2tra import SpecialNoMateDiffBNDPairIndependentToTRAConverter
from octopusv.converter.snmd_dndpr_tra2tra import SpecialNoMateDiffBNDPairReciprocalTranslocationToTRAConverter
from octopusv.converter.stra2tra import SingleTRAToTRAConverter
from octopusv.transformer.mp_bnd import MatePairBNDTransformer
from octopusv.transformer.no_bnd import NonBNDTransformer
from octopusv.transformer.same_chr_sv import SameChrSVTransformer
from octopusv.transformer.snmd_bndp import SpecialNoMateDiffBNDPairTransformer
from octopusv.transformer.stra import SingleTRATransformer
from octopusv.utils.normal_vcf_parser import parse_vcf
from octopusv.utils.svcf_utils import write_sv_vcf


def correct(
        input_vcf: Path | None = typer.Argument(
            None, exists=True, dir_okay=False, resolve_path=True, help="Input VCF file to correct."
        ),
        output: Path | None = typer.Argument(None, dir_okay=False, resolve_path=True, help="Output file path."),
        input_option: Path | None = typer.Option(
            None, "--input-file", "-i", exists=True, dir_okay=False, resolve_path=True,
            help="Input VCF file to correct."
        ),
        output_option: Path | None = typer.Option(
            None, "--output-file", "-o", dir_okay=False, resolve_path=True, help="Output file path."
        ),
        pos_tolerance: int = typer.Option(
            3,
            "--pos-tolerance",
            "-pt",
            help="Position tolerance for identifying mate BND events, default=3, recommend not to set larger than 5",
        ),
        # Quality filtering parameters
        min_qual: float | None = typer.Option(
            None, "--min-qual", help="Minimum QUAL score to keep variants"
        ),
        max_qual: float | None = typer.Option(
            None, "--max-qual", help="Maximum QUAL score to keep variants"
        ),
        min_support: int | None = typer.Option(
            None, "--min-support", help="Minimum supporting reads to keep variants"
        ),
        max_support: int | None = typer.Option(
            None, "--max-support", help="Maximum supporting reads to keep variants"
        ),
        min_depth: int | None = typer.Option(
            None, "--min-depth", help="Minimum total depth to keep variants"
        ),
        max_depth: int | None = typer.Option(
            None, "--max-depth", help="Maximum total depth to keep variants"
        ),
        min_gq: int | None = typer.Option(
            None, "--min-gq", help="Minimum genotype quality to keep variants"
        ),
        min_svlen: int | None = typer.Option(
            None, "--min-svlen", help="Minimum SV length to keep variants"
        ),
        max_svlen: int | None = typer.Option(
            None, "--max-svlen", help="Maximum SV length to keep variants"
        ),
        filter_pass: bool = typer.Option(
            False, "--filter-pass", help="Only keep variants with FILTER=PASS"
        ),
        exclude_nocall: bool = typer.Option(
            False, "--exclude-nocall", help="Exclude variants with ./. genotype"
        ),
):
    """Correct SV events with optional quality filtering."""

    # Determine input file
    if input_vcf and input_option:
        typer.echo(
            "Error: Please specify input file either as an argument or with -i/--input-file, not both.", err=True
        )
        raise typer.Exit(code=1)
    if input_vcf:
        input_file = input_vcf
    elif input_option:
        input_file = input_option
    else:
        typer.echo("Error: Input file is required.", err=True)
        raise typer.Exit(code=1)

    # Determine output file
    if output and output_option:
        typer.echo(
            "Error: Please specify output file either as an argument or with -o/--output-file, not both.", err=True
        )
        raise typer.Exit(code=1)
    if output:
        output_file = output
    elif output_option:
        output_file = output_option
    else:
        typer.echo("Error: Output file is required.", err=True)
        raise typer.Exit(code=1)

    # Parse the input VCF file
    # non_bnd_events means DEL, INV, INS, DUP
    contig_lines, same_chr_bnd_events, diff_chr_bnd_events, non_bnd_events = parse_vcf(input_file)

    # Initialize quality filter if any filtering parameters are provided
    filter_params = [min_qual, max_qual, min_support, max_support, min_depth, max_depth,
                     min_gq, min_svlen, max_svlen, filter_pass, exclude_nocall]
    has_quality_filters = any([
        min_qual is not None, max_qual is not None,
        min_support is not None, max_support is not None,
        min_depth is not None, max_depth is not None,
        min_gq is not None, min_svlen is not None, max_svlen is not None,
        filter_pass, exclude_nocall
    ])

    if has_quality_filters:
        # Import QualityFilter only when needed
        from octopusv.filter import QualityFilter

        quality_filter = QualityFilter(
            min_qual=min_qual,
            max_qual=max_qual,
            min_support=min_support,
            max_support=max_support,
            min_depth=min_depth,
            max_depth=max_depth,
            min_gq=min_gq,
            min_svlen=min_svlen,
            max_svlen=max_svlen,
            filter_pass=filter_pass,
            exclude_nocall=exclude_nocall
        )

        # Apply quality filtering to all event types
        def apply_quality_filter(events):
            """Apply quality filter to a list of events."""
            return [event for event in events if quality_filter.filter_event(event)]

        # Filter all event categories
        same_chr_bnd_events = apply_quality_filter(same_chr_bnd_events)
        diff_chr_bnd_events = apply_quality_filter(diff_chr_bnd_events)
        non_bnd_events = apply_quality_filter(non_bnd_events)

        # Print filtering statistics
        quality_filter.print_stats()
        typer.echo(f"Quality filtering completed.")

    # Extract mate BND and no mate events, they are all with different chromosomes
    mate_bnd_pairs = find_mate_bnd_events(diff_chr_bnd_events, pos_tolerance=pos_tolerance)
    no_mate_events = find_no_mate_events(diff_chr_bnd_events, pos_tolerance=pos_tolerance)

    # Further classify no_mate_events into special_no_mate_diff_bnd_pair and other_single_TRA
    special_no_mate_diff_bnd_pairs, other_single_TRA = find_special_no_mate_diff_bnd_pair_and_other_single_tra(
        no_mate_events,
        pos_tolerance=pos_tolerance,
    )

    # Initialize the EventTransformer with a list of transform strategies for each type of events
    same_chr_sv_transformer = SameChrSVTransformer(
        [
            BNDKeepingConverter(),
        ]
    )

    mate_bnd_pair_transformer = MatePairBNDTransformer(
        [
            MatePairReciprocalTranslocationToTRAConverter(),
            MatePairIndependentToTRAConverter(),
            MatePairMergeToTRAConverter(),
        ]
    )

    special_no_mate_diff_bnd_pair_transformer = SpecialNoMateDiffBNDPairTransformer(
        [
            SpecialNoMateDiffBNDPairReciprocalTranslocationToTRAConverter(),
            SpecialNoMateDiffBNDPairIndependentToTRAConverter(),
        ]
    )

    single_TRA_transformer = SingleTRATransformer([SingleTRAToTRAConverter()])
    non_bnd_transformer = NonBNDTransformer([NonBNDConverter()])

    # Apply all transformation strategies to the events
    same_chr_sv_transformed_events = same_chr_sv_transformer.apply_transforms(same_chr_bnd_events)
    mate_pair_transformed_events = mate_bnd_pair_transformer.apply_transforms(mate_bnd_pairs)
    special_no_mate_diff_bnd_pair_transformed_events = special_no_mate_diff_bnd_pair_transformer.apply_transforms(
        special_no_mate_diff_bnd_pairs
    )
    single_TRA_transformed_events = single_TRA_transformer.apply_transforms(other_single_TRA)
    non_bnd_transformed_events = non_bnd_transformer.apply_transforms(non_bnd_events)

    # Merge all transformed events
    all_transformed_events = (
            same_chr_sv_transformed_events
            + mate_pair_transformed_events
            + special_no_mate_diff_bnd_pair_transformed_events
            + single_TRA_transformed_events
            + non_bnd_transformed_events
    )

    # Write the transformed events to the output file
    write_sv_vcf(contig_lines, all_transformed_events, output_file)
    typer.echo(f"Corrected SV events written to {output_file}")


def find_mate_bnd_events(events, pos_tolerance=3):
    """Find mate BND events."""
    event_dict = {}
    mate_bnd_pairs = []

    for event in events:
        try:
            chrom_alt, pos_alt = get_alt_chrom_pos(event.alt)
            if chrom_alt is None or pos_alt is None:
                continue

            key_for_searching = (event.chrom, event.pos, chrom_alt, pos_alt)

            # Generate all possible reverse keys
            possible_reverse_keys = []
            for i in range(-pos_tolerance, pos_tolerance + 1):
                for j in range(-pos_tolerance, pos_tolerance + 1):
                    try:
                        new_key = (chrom_alt, pos_alt + i, event.chrom, event.pos + j)
                        possible_reverse_keys.append(new_key)
                    except TypeError:
                        continue

            for reverse_key in possible_reverse_keys:
                if reverse_key in event_dict:
                    mate_bnd_pairs.append((event_dict.pop(reverse_key), event))
                    break
            else:
                event_dict[key_for_searching] = event

        except (TypeError, ValueError) as e:
            logging.info(f"Error processing event: {e!s}")
            continue

    return mate_bnd_pairs


def find_no_mate_events(events, pos_tolerance=3):
    """Find no mate events."""
    event_dict = {}
    mate_bnd_pairs = []

    for event in events:
        try:
            chrom_alt, pos_alt = get_alt_chrom_pos(event.alt)
            if chrom_alt is None or pos_alt is None:
                continue

            key = (event.chrom, event.pos, chrom_alt, pos_alt)

            # Generate all possible reverse keys
            possible_reverse_keys = []
            for i in range(-pos_tolerance, pos_tolerance + 1):
                for j in range(-pos_tolerance, pos_tolerance + 1):
                    try:
                        new_key = (chrom_alt, pos_alt + i, event.chrom, event.pos + j)
                        possible_reverse_keys.append(new_key)
                    except TypeError:
                        continue

            mate_found = False
            for reverse_key in possible_reverse_keys:
                if reverse_key in event_dict:
                    for mate_event in event_dict[reverse_key]:
                        mate_bnd_pairs.append((mate_event, event))
                        event_dict[reverse_key].remove(mate_event)
                        mate_found = True
                        break
                    if mate_found:
                        if len(event_dict[reverse_key]) == 0:
                            del event_dict[reverse_key]
                        break

            if not mate_found:
                if key not in event_dict:
                    event_dict[key] = []
                event_dict[key].append(event)

        except (TypeError, ValueError) as e:
            logging.info(f"Error processing event: {e!s}")
            continue

    return [single_event for event_list in event_dict.values() for single_event in event_list]


def find_special_no_mate_diff_bnd_pair_and_other_single_tra(events, pos_tolerance=3):
    """Extract special no mate diff bnd pair and other single TRA events."""
    event_dict = {}
    special_no_mate_diff_bnd_pair = []
    other_single_TRA = []

    for event in events:
        try:
            chrom_alt, pos_alt = get_alt_chrom_pos(event.alt)
            if chrom_alt is None or pos_alt is None:
                continue

            key = (event.chrom, event.pos, chrom_alt, pos_alt)

            # Generate all possible reverse keys
            possible_reverse_keys = []
            for i in range(-pos_tolerance, pos_tolerance + 1):
                for j in range(-pos_tolerance, pos_tolerance + 1):
                    try:
                        new_key = (event.chrom, event.pos + i, chrom_alt, pos_alt + j)
                        possible_reverse_keys.append(new_key)
                    except TypeError:
                        continue

            mate_found = False
            for reverse_key in possible_reverse_keys:
                if reverse_key in event_dict:
                    special_no_mate_diff_bnd_pair.append((event_dict.pop(reverse_key), event))
                    mate_found = True
                    break

            if not mate_found:
                event_dict[key] = event

        except (TypeError, ValueError) as e:
            logging.info(f"Error processing event: {e!s}")
            continue

    other_single_TRA = list(event_dict.values())
    return special_no_mate_diff_bnd_pair, other_single_TRA