import os
import sys
from pathlib import Path

import typer

from octopusv.merger.name_mapper import NameMapper
from octopusv.merger.sv_merger import SVMerger
from octopusv.merger.upset_plotter import UpSetPlotter
from octopusv.utils.SV_classifier_by_chromosome import SVClassifiedByChromosome
from octopusv.utils.SV_classifier_by_type import SVClassifierByType
from octopusv.utils.svcf_parser import SVCFFileEventCreator


def _echo(message: str = "") -> None:
    """Print human-facing CLI messages to stderr.

    This keeps stdout clean for future machine-readable output or pipelines.
    """
    typer.echo(message, err=True)


def _default_labels_from_input_files(input_files: list[Path | str]) -> list[str]:
    """Return default labels inferred from input file basenames."""
    return [os.path.splitext(os.path.basename(str(file)))[0] for file in input_files]


def _format_labels(labels: list[str]) -> str:
    """Format labels for human-readable CLI output."""
    return ", ".join(labels) if labels else "."


def _format_input_order(input_files: list[Path | str]) -> str:
    """Format input file order for human-readable CLI output."""
    return str([str(file) for file in input_files])


def _labels_match_default(labels: list[str], input_files: list[Path | str]) -> bool:
    """Return whether user-provided labels match default basename labels."""
    return labels == _default_labels_from_input_files(input_files)


def _normalize_input_path(path: Path | str) -> str:
    """Return a stable real-path identity for one merge input."""
    return os.path.normcase(
        os.path.realpath(
            os.path.abspath(str(path))
        )
    )


def _validate_distinct_input_files(input_files: list[Path | str]) -> None:
    """Reject two merge arguments that resolve to the same physical file."""
    seen = {}

    for input_file in input_files:
        normalized = _normalize_input_path(input_file)
        previous = seen.get(normalized)
        if previous is not None:
            raise ValueError(
                "Two merge inputs resolve to the same physical file: "
                f"{previous!r} and {str(input_file)!r}. "
                "Each merge input must refer to a distinct file."
            )
        seen[normalized] = str(input_file)


def _validate_unique_display_labels(
    *,
    labels: list[str],
    input_files: list[Path | str],
    mode: str,
) -> None:
    """Reject display labels that would make two inputs indistinguishable."""
    label_to_files = {}
    for label, input_file in zip(labels, input_files):
        label_to_files.setdefault(label, []).append(str(input_file))

    duplicates = {
        label: files
        for label, files in label_to_files.items()
        if len(files) > 1
    }
    if not duplicates:
        return

    option = "--sample-names" if mode == "sample" else "--caller-names"
    details = "; ".join(
        f"{label!r}: {files}"
        for label, files in duplicates.items()
    )
    raise ValueError(
        "Merge input labels must be unique. "
        f"Duplicate label(s): {details}. "
        f"Use {option} to provide one unique label per input file."
    )


def _preflight_svcf_shape(input_file: Path | str, mode: str) -> None:
    """Reject input layouts that the selected merge mode cannot preserve.

    This is intentionally a lightweight byte-level scan: it checks only the
    #CHROM layout and tab counts, without splitting or parsing full records.

    Caller mode requires every input record to carry exactly one trailing
    evidence/sample column. Re-merging a variable-width caller-merged SVCF
    would otherwise discard nested evidence.

    Sample mode intentionally allows a caller-merged SVCF with one header
    sample label and multiple caller evidence blocks per record. This is the
    supported per-sample-caller-merge -> population-merge workflow. It does
    not allow a true multi-sample SVCF to be nested as one input file.
    """
    found_chrom_header = False

    with open(input_file, "rb") as handle:
        for line_number, line in enumerate(handle, 1):
            if line.startswith(b"#CHROM"):
                found_chrom_header = True
                # A valid SVCF has nine fixed columns through FORMAT.
                # fields = tabs + 1, so trailing columns = tabs - 8.
                header_sample_count = max(0, line.count(b"\t") - 8)

                if mode == "sample" and header_sample_count > 1:
                    raise ValueError(
                        f"Input {str(input_file)!r} has {header_sample_count} "
                        "sample columns. `merge --mode sample` currently "
                        "expects one input file per biological sample. "
                        "Split or subset the file to one sample per input "
                        "before merging."
                    )
                continue

            if line.startswith(b"#") or not line.strip():
                continue

            if not found_chrom_header:
                raise ValueError(
                    f"Input {str(input_file)!r} is missing a #CHROM header "
                    "before its data records."
                )

            if mode == "caller":
                evidence_count = max(0, line.count(b"\t") - 8)
                if evidence_count > 1:
                    raise ValueError(
                        f"Input {str(input_file)!r} contains {evidence_count} "
                        f"evidence/sample columns on data line {line_number}. "
                        "Caller-mode re-merge of multi-evidence SVCF records "
                        "is not supported because nested evidence cannot be "
                        "preserved unambiguously. Use single-evidence SVCF "
                        "inputs for `merge --mode caller`."
                    )
            else:
                # In sample mode the header determines whether this is a true
                # multi-sample input. Variable-width caller evidence below a
                # single SAMPLE header is intentionally allowed.
                break

    if not found_chrom_header:
        raise ValueError(
            f"Input {str(input_file)!r} is missing the required #CHROM header."
        )


def _preflight_merge_inputs(
    *,
    input_files: list[Path | str],
    labels: list[str],
    mode: str,
) -> None:
    """Run merge-specific identity and input-shape checks before parsing."""
    _validate_distinct_input_files(input_files)
    _validate_unique_display_labels(
        labels=labels,
        input_files=input_files,
        mode=mode,
    )

    for input_file in input_files:
        _preflight_svcf_shape(input_file, mode)


def _mark_sample_input_collapses(events) -> None:
    """Preserve sample-input evidence and annotate current collapse counts.

    A caller-merged SVCF may legitimately be used as one biological-sample
    input to a population/sample-mode merge. Until the 1.0 sample-consensus
    synthesis is connected, the existing writer still uses the first evidence
    block as the sample-level representation.

    This function therefore has two intentionally separate responsibilities:

    1. Preserve the complete caller-evidence payload on the existing sample
       dictionary under a private ``_octopusv_`` key so it survives selection
       and reaches the sample-mode writer without changing merge semantics.
    2. Keep the historical private collapse counter unchanged so current
       output/warnings remain byte-for-byte compatible during this transport
       step.

    Runtime-only ``_octopusv_`` keys are not part of FORMAT and must never be
    serialized into SVCF output.
    """
    collapse_key = "_octopusv_collapsed_evidence_count"
    payload_key = "_octopusv_evidence_payload"

    for event in events:
        raw_columns = list(getattr(event, "raw_sample_columns", []) or [])
        sample_data = getattr(event, "sample", None)
        if not isinstance(sample_data, dict):
            continue

        info = getattr(event, "info", None)
        if not isinstance(info, dict):
            info = {}

        # Transport only: keep the raw evidence and its explicit positional
        # source fields intact. Do not interpret, collapse, or vote here.
        sample_data[payload_key] = {
            "format": getattr(event, "format", ""),
            "blocks": tuple(raw_columns),
            "sources": info.get("SOURCES"),
            "source_ids": info.get("SOURCE_IDS"),
        }

        additional = max(0, len(raw_columns) - 1)
        if additional == 0:
            continue

        existing = sample_data.get(collapse_key, 0)
        try:
            existing = int(existing)
        except (TypeError, ValueError):
            existing = 0

        sample_data[collapse_key] = max(0, existing) + additional


def _describe_merge_rule(
    *,
    expression: str | None,
    intersect: bool,
    union: bool,
    specific: list[Path] | None,
    min_support: int | None,
    exact_support: int | None,
    max_support: int | None,
) -> list[str]:
    """Return a short explanation of the selected merge strategy."""
    lines = []

    if expression:
        lines.append(f"Merge rule: --expression {expression}")
        lines.append("Only records satisfying the logical source expression were retained.")
        return lines

    if intersect:
        lines.append("Merge rule: --intersect")
        lines.append("Only records supported by all input files were retained.")
        return lines

    if union:
        lines.append("Merge rule: --union")
        lines.append("Records supported by at least one input file were retained.")
        return lines

    if specific:
        specific_files = _format_input_order(specific)
        lines.append(f"Merge rule: --specific {specific_files}")
        lines.append("Only records specific to the selected input file(s) were retained.")
        return lines

    if exact_support is not None:
        lines.append(f"Merge rule: --exact-support {exact_support}")
        lines.append(f"Only records supported by exactly {exact_support} input file(s) were retained.")
        return lines

    if min_support is not None and max_support is not None:
        lines.append(f"Merge rule: --min-support {min_support} --max-support {max_support}")
        lines.append(
            f"Records supported by fewer than {min_support} or more than {max_support} "
            "input file(s) were excluded."
        )
        return lines

    if min_support is not None:
        lines.append(f"Merge rule: --min-support {min_support}")
        lines.append(f"Records supported by fewer than {min_support} input file(s) were excluded.")
        return lines

    if max_support is not None:
        lines.append(f"Merge rule: --max-support {max_support}")
        lines.append(f"Records supported by more than {max_support} input file(s) were excluded.")
        return lines

    return lines


def _print_initial_configuration(
    *,
    mode: str,
    labels: list[str],
    input_files: list[Path | str],
) -> None:
    """Print a concise configuration summary before running merge."""
    if mode == "sample":
        _echo("Info: Running in sample mode.")
        _echo(f"Info: Sample labels: {_format_labels(labels)}")
        _echo(f"Info: Input file order: {_format_input_order(input_files)}")
    else:
        _echo("Info: Running in caller mode.")
        _echo(f"Info: Caller labels: {_format_labels(labels)}")
        _echo(f"Info: Input file order: {_format_input_order(input_files)}")


def _print_name_mapping_note(
    *,
    mode: str,
    labels: list[str],
    input_files: list[Path | str],
    caller_names: str | None,
    sample_names: str | None,
) -> None:
    """Explain what --caller-names or --sample-names does and does not do."""
    default_labels = _default_labels_from_input_files(input_files)

    if mode == "caller":
        if caller_names:
            if labels == default_labels:
                _echo(
                    "Note: Provided --caller-names match the default input basenames, "
                    "so output labels are unchanged."
                )
            else:
                _echo(
                    "Note: --caller-names changes caller labels in SOURCES and downstream "
                    "OctopuSV selection; it does not change which SVs are merged."
                )
        else:
            _echo("Note: No --caller-names provided; using input file basenames as caller labels.")
        return

    if mode == "sample":
        if sample_names:
            if labels == default_labels:
                _echo(
                    "Note: Provided --sample-names match the default input basenames, "
                    "so output sample column labels are unchanged."
                )
            else:
                _echo(
                    "Note: --sample-names changes output sample column labels; "
                    "it does not change which SVs are merged."
                )
        else:
            _echo("Note: No --sample-names provided; using input file basenames as sample labels.")


def _print_success_message(
    *,
    output_file: Path,
    result_count: int,
    input_count: int,
    mode: str,
    labels: list[str],
    input_files: list[Path | str],
    caller_names: str | None,
    sample_names: str | None,
    expression: str | None,
    intersect: bool,
    union: bool,
    specific: list[Path] | None,
    min_support: int | None,
    exact_support: int | None,
    max_support: int | None,
) -> None:
    """Print a concise, mode-aware merge summary."""
    _echo(f"Successfully merged {result_count} events from {input_count} input files.")
    _echo(f"Merged results written to {output_file}")
    _echo("")

    if mode == "caller":
        _echo("Mode: caller")
        _echo("Output layout: one SAMPLE column with multiple caller evidence blocks.")
        _echo(f"Caller labels: {_format_labels(labels)}")
        _echo("SOURCES, SOURCE_IDS, and evidence blocks follow the input file order.")
    else:
        _echo("Mode: sample")
        _echo(f"Output layout: {len(labels)} sample/input columns: {_format_labels(labels)}")
        _echo("SOURCES and SOURCE_IDS follow the retained sample/input order.")

    _print_name_mapping_note(
        mode=mode,
        labels=labels,
        input_files=input_files,
        caller_names=caller_names,
        sample_names=sample_names,
    )

    merge_rule_lines = _describe_merge_rule(
        expression=expression,
        intersect=intersect,
        union=union,
        specific=specific,
        min_support=min_support,
        exact_support=exact_support,
        max_support=max_support,
    )
    if merge_rule_lines:
        _echo("")
        for line in merge_rule_lines:
            _echo(line)

    _echo("")
    _echo(
        "Output type: SVCF. Convert with `octopusv svcf2vcf` before using "
        "standard VCF tools such as bcftools/vcftools."
    )


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
                    line = line.strip()
                    if line.startswith("##contig=<") and line.endswith(">"):
                        content = line[len("##contig=<"): -1]
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

        # Mode parameters
        mode: str = typer.Option(
            "caller",
            "--mode",
            help=(
                "Merge mode: 'caller' for the same sample analyzed by different callers; "
                "'sample' for different samples or input datasets."
            ),
        ),
        caller_names: str = typer.Option(
            None,
            "--caller-names",
            help=(
                "Comma-separated display labels for input callers in caller mode. "
                "Defaults to input file basenames. This changes SOURCES labels, "
                "not merge criteria."
            ),
        ),
        sample_names: str = typer.Option(
            None,
            "--sample-names",
            help=(
                "Comma-separated display labels for input samples in sample mode. "
                "Defaults to input file basenames. This changes output sample column names, "
                "not merge criteria."
            ),
        ),

        # Existing merge strategy parameters
        intersect: bool = typer.Option(False, "--intersect", help="Apply intersection strategy for merging."),
        union: bool = typer.Option(False, "--union", help="Apply union strategy for merging."),
        specific: list[Path] = typer.Option(
            None, "--specific", help="Extract SVs that are specifically supported by provided files."
        ),
        min_support: int = typer.Option(None, "--min-support", help="Minimum number of files that must support an SV."),
        exact_support: int = typer.Option(
            None, "--exact-support", help="Exact number of files that must support an SV."
        ),
        max_support: int = typer.Option(None, "--max-support", help="Maximum number of files that can support an SV."),
        expression: str = typer.Option(
            None,
            "--expression",
            help="Logical expression for complex file combinations (e.g., '(A AND B) AND NOT (C OR D)')",
        ),

        # Ordinary SV matching parameters
        max_distance: int | None = typer.Option(
            None,
            "--max-distance",
            help=(
                "Override the maximum breakpoint distance used to merge ordinary SVs. "
                "If not provided, OctopuSV uses SV-type- and size-aware adaptive thresholds."
            ),
        ),
        max_length_ratio: float | None = typer.Option(
            None,
            "--max-length-ratio",
            help=(
                "Override the maximum SV length ratio used to merge ordinary SVs. "
                "If not provided, OctopuSV uses SV-type-specific thresholds."
            ),
        ),
        min_jaccard: float = typer.Option(
            0.0,
            "--min-jaccard",
            help=(
                "Minimum interval Jaccard overlap required for DEL, DUP, and INV merging. "
                "Default: 0 (disabled). "
                "INS, TRA, and BND are not evaluated with interval Jaccard."
            ),
        ),
        tra_delta: int = typer.Option(
            50, "--tra-delta", help="Position uncertainty threshold for TRA events (in base pairs)."
        ),
        tra_min_overlap_ratio: float = typer.Option(
            0.5, "--tra-min-overlap", help="Minimum overlap ratio for TRA events."
        ),
        tra_strand_consistency: bool = typer.Option(
            True, "--tra-strand-consistency", help="Whether to require strand consistency for TRA events."
        ),
        bnd_delta: int = typer.Option(
            50, "--bnd-delta", help="Position uncertainty threshold for BND events (in base pairs)."
        ),

        # Visualization parameters
        upsetr: bool = typer.Option(
            False, "--upsetr", help="Generate UpSet plot visualization of input file intersections."
        ),
        upsetr_output: Path = typer.Option(
            None,
            "--upsetr-output",
            help="Output path for UpSet plot. If not provided, will use output_file basename with _upset.png suffix.",
        ),
):
    """Merge multiple SVCF files based on the selected strategy."""

    # Validate mode parameter.
    if mode not in ["caller", "sample"]:
        typer.echo("Error: --mode must be either 'caller' or 'sample'.", err=True)
        raise typer.Exit(code=1)

    # Validate mode-specific parameters.
    if mode == "caller" and sample_names:
        typer.echo("Error: --sample-names can only be used with --mode sample.", err=True)
        raise typer.Exit(code=1)

    if mode == "sample" and caller_names:
        typer.echo("Error: --caller-names can only be used with --mode caller.", err=True)
        raise typer.Exit(code=1)

    # Handle input files while preserving user-provided order.
    all_input_files = []

    using_i_option = "-i" in sys.argv or "--input-file" in sys.argv

    if using_i_option:
        if input_option:
            all_input_files.extend(input_option)
        if input_files:
            all_input_files.extend(input_files)
    else:
        if input_files:
            all_input_files.extend(input_files)
        if input_option:
            all_input_files.extend(input_option)

    if not all_input_files:
        typer.echo("Error: No input files provided.", err=True)
        raise typer.Exit(code=1)

    if specific and not specific[0]:
        typer.echo("Error: --specific option requires at least one file.", err=True)
        raise typer.Exit(code=1)

    if min_support is not None and min_support < 1:
        typer.echo("Error: --min-support must be a positive integer.", err=True)
        raise typer.Exit(code=1)

    # Validate ordinary SV matching parameters.
    if max_distance is not None and max_distance < 0:
        typer.echo("Error: --max-distance must be non-negative.", err=True)
        raise typer.Exit(code=1)

    if max_length_ratio is not None and max_length_ratio < 1.0:
        typer.echo("Error: --max-length-ratio must be at least 1.0.", err=True)
        raise typer.Exit(code=1)

    if not 0.0 <= min_jaccard <= 1.0:
        typer.echo("Error: --min-jaccard must be between 0 and 1.", err=True)
        raise typer.Exit(code=1)

    # Build name mapper.
    name_mapper = None
    try:
        if mode == "caller" and caller_names:
            names = [name.strip() for name in caller_names.split(",")]
            if len(names) != len(all_input_files):
                typer.echo(
                    f"Error: Number of caller names ({len(names)}) must match "
                    f"number of input files ({len(all_input_files)}).",
                    err=True,
                )
                raise typer.Exit(code=1)
            name_mapper = NameMapper(all_input_files, mode="caller", custom_names=names)

        elif mode == "sample":
            names = None
            if sample_names:
                names = [name.strip() for name in sample_names.split(",")]
                if len(names) != len(all_input_files):
                    typer.echo(
                        f"Error: Number of sample names ({len(names)}) must match "
                        f"number of input files ({len(all_input_files)}).",
                        err=True,
                    )
                    raise typer.Exit(code=1)
            name_mapper = NameMapper(all_input_files, mode="sample", custom_names=names)

        elif mode == "caller":
            name_mapper = NameMapper(all_input_files, mode="caller")

    except ValueError as e:
        typer.echo(f"Error: {e}", err=True)
        raise typer.Exit(code=1)

    labels = name_mapper.get_all_display_names() if name_mapper else []

    try:
        _preflight_merge_inputs(
            input_files=all_input_files,
            labels=labels,
            mode=mode,
        )
    except (OSError, ValueError) as e:
        _echo(f"Error: {e}")
        raise typer.Exit(code=1)

    _print_initial_configuration(mode=mode, labels=labels, input_files=all_input_files)

    # Get contig information from input files.
    input_filenames = [str(file) for file in all_input_files]
    contigs = get_contigs_from_svcf(input_filenames)

    # Process SV events.
    sv_event_creator = SVCFFileEventCreator(input_filenames)
    sv_event_creator.parse()

    if mode == "sample":
        _mark_sample_input_collapses(sv_event_creator.events)

    classifier = SVClassifierByType(sv_event_creator.events)
    classifier.classify()

    chromosome_classifier = SVClassifiedByChromosome(classifier.get_classified_events())
    chromosome_classifier.classify()

    # Initialize merger with preserved file order.
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

    # Apply merge strategy.
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

    # Write merged results.
    write_summary = sv_merger.write_results(
        output_file,
        results,
        contigs,
        mode,
        name_mapper,
        input_filenames,
    ) or {}

    if mode == "sample":
        collapse_records = int(write_summary.get("sample_collapse_records", 0) or 0)
        collapse_blocks = int(write_summary.get("sample_collapsed_evidence_blocks", 0) or 0)
        if collapse_records > 0:
            _echo(
                "Warning: Sample-mode output collapsed multiple evidence blocks "
                f"to one deterministic block per sample in {collapse_records} "
                f"merged record(s); {collapse_blocks} additional evidence block(s) "
                "were not serialized as separate sample columns. The retained "
                "block is the first deterministic evidence block from each input; "
                "for caller-merged per-sample inputs, this is the first evidence "
                "column. If caller input order differs across samples, the retained "
                "caller may also differ. The merged events were retained."
            )

    _print_success_message(
        output_file=output_file,
        result_count=len(results),
        input_count=len(all_input_files),
        mode=mode,
        labels=labels,
        input_files=all_input_files,
        caller_names=caller_names,
        sample_names=sample_names,
        expression=expression,
        intersect=intersect,
        union=union,
        specific=specific,
        min_support=min_support,
        exact_support=exact_support,
        max_support=max_support,
    )

    # Generate UpSet plot if requested.
    if upsetr:
        try:
            plot_file = str(upsetr_output) if upsetr_output else str(output_file).rsplit(".", 1)[0] + "_upset.png"
            plotter = UpSetPlotter(sv_merger.get_all_merged_events(), all_input_files)
            plotter.plot(plot_file)
            _echo(f"UpSet plot written to {plot_file}")
        except ImportError:
            typer.echo(
                "Warning: Could not generate UpSet plot. Please ensure matplotlib and numpy are installed.",
                err=True,
            )
        except Exception as e:
            typer.echo(f"Warning: Failed to generate UpSet plot: {e!s}", err=True)


if __name__ == "__main__":
    typer.run(merge)
