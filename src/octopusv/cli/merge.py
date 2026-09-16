import os
import sys
from pathlib import Path

import typer

from octopusv.merger.name_mapper import NameMapper, default_label_from_path
from octopusv.merger.sv_merge_selection import validate_selection_inputs
from octopusv.merger.sv_merger import SVMerger
from octopusv.merger.upset_plotter import UpSetPlotter
from octopusv.utils.SV_classifier_by_chromosome import SVClassifiedByChromosome
from octopusv.utils.SV_classifier_by_type import SVClassifierByType
from octopusv.utils.atomic_write import atomic_output_path
from octopusv.utils.source_path import normalize_source_path
from octopusv.utils.text_io import open_text_auto
from octopusv.utils.svcf_parser import SVCFFileEventCreator
from octopusv.utils.svcf_schema import (
    MODE_CALLER,
    MODE_MULTI,
    parse_identity_from_meta_lines,
    validate_format_for_identity,
    validate_header_columns_for_identity,
    validate_source_label,
    validate_versioned_identity,
)


def _echo(message: str = "") -> None:
    """Print human-facing CLI messages to stderr.

    This keeps stdout clean for future machine-readable output or pipelines.
    """
    typer.echo(message, err=True)


def _default_labels_from_input_files(input_files: list[Path | str]) -> list[str]:
    """Return default labels inferred from input file basenames."""
    return [default_label_from_path(str(file)) for file in input_files]


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
    return normalize_source_path(path)



def _read_merge_input_list(list_path: Path | str) -> tuple[list[Path], list[str] | None]:
    """Read a cohort input list.

    Each non-comment line is either ``PATH`` or ``PATH<TAB>LABEL``. Relative
    paths are resolved relative to the list file, not the current working
    directory. A labeled list must label every data line.
    """
    list_path = Path(os.path.abspath(os.path.expanduser(str(list_path))))
    base_dir = list_path.parent
    paths: list[Path] = []
    labels: list[str | None] = []

    with open_text_auto(list_path) as handle:
        for line_number, raw in enumerate(handle, 1):
            stripped = raw.strip()
            if not stripped or stripped.startswith("#"):
                continue

            parts = raw.rstrip("\r\n").split("\t")
            if len(parts) > 2:
                raise ValueError(
                    f"Invalid --input-list line {line_number}: expected PATH or "
                    "PATH<TAB>LABEL."
                )

            raw_path = parts[0].strip()
            if not raw_path:
                raise ValueError(
                    f"Invalid --input-list line {line_number}: path is empty."
                )

            candidate = Path(raw_path).expanduser()
            if not candidate.is_absolute():
                candidate = base_dir / candidate
            # Preserve the user-visible path (including a symlink basename) for
            # caller/sample identity. Physical-file duplicate detection uses
            # normalize_source_path/realpath separately.
            candidate = Path(os.path.abspath(candidate))

            if not candidate.exists():
                raise ValueError(
                    f"Input file listed on line {line_number} does not exist: "
                    f"{candidate}"
                )

            label = None
            if len(parts) == 2:
                label = parts[1].strip()
                if not label:
                    raise ValueError(
                        f"Invalid --input-list line {line_number}: label is empty."
                    )

            paths.append(candidate)
            labels.append(label)

    if not paths:
        raise ValueError(f"Input list {str(list_path)!r} contains no input files.")

    has_labels = [label is not None for label in labels]
    if any(has_labels) and not all(has_labels):
        raise ValueError(
            "A labeled --input-list must provide a TAB-separated label for "
            "every input row."
        )

    return paths, [str(label) for label in labels] if all(has_labels) else None


def _validate_display_labels(
    *,
    labels: list[str],
    mode: str,
) -> None:
    """Validate final caller/sample labels before any expensive merge work."""
    option = "--sample-names" if mode == "sample" else "--caller-names"
    for label in labels:
        try:
            validate_source_label(label)
        except ValueError as exc:
            raise ValueError(
                f"Invalid {mode} label {label!r}: {exc} "
                f"Use {option} or a labeled --input-list to provide a valid label."
            ) from exc

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


def _extract_contig_lengths_from_meta_lines(
    meta_lines: list[str],
    *,
    input_file: Path | str,
) -> dict[str, str]:
    """Return explicit numeric ##contig lengths declared by one input."""
    contigs: dict[str, str] = {}
    for line in meta_lines:
        if not (line.startswith("##contig=<") and line.endswith(">")):
            continue

        content = line[len("##contig=<") : -1]
        contig_id = None
        contig_length = None
        for part in content.split(","):
            if part.startswith("ID="):
                contig_id = part.split("=", 1)[1]
            elif part.startswith("length="):
                raw_length = part.split("=", 1)[1]
                try:
                    int(raw_length)
                except ValueError:
                    # A non-numeric declaration is not a usable known length for
                    # conflict detection. The validator remains responsible for
                    # broader header/schema diagnostics.
                    raw_length = None
                contig_length = raw_length

        if not contig_id or contig_length is None:
            continue

        previous = contigs.get(contig_id)
        if previous is not None and int(previous) != int(contig_length):
            raise ValueError(
                f"Input {str(input_file)!r} declares conflicting lengths for "
                f"contig {contig_id!r}: {previous} and {contig_length}."
            )
        contigs[contig_id] = contig_length

    return contigs


def _preflight_svcf_shape(input_file: Path | str, mode: str) -> dict[str, str]:
    """Reject merge inputs that are not caller-evidence SVCF.

    Merge consumes caller evidence. It may accept legacy unversioned SVCF or
    explicit SVCF 1.1 caller-mode files, but it must never accept synthesized
    sample/multi SVCF as caller evidence. FORMAT schema is checked directly so
    stripping a mode marker cannot make synthesized sample calls mergeable.

    Returns explicit numeric contig lengths found in the header so cross-input
    reference conflicts can be rejected before parsing/merging begins.
    """
    found_chrom_header = False
    meta_lines: list[str] = []
    identity = None
    header_sample_count = None
    contig_lengths: dict[str, str] = {}

    with open_text_auto(input_file) as handle:
        for line_number, line in enumerate(handle, 1):
            if line.startswith("##"):
                meta_lines.append(line.rstrip("\r\n"))
                continue

            if line.startswith("#CHROM"):
                try:
                    identity = parse_identity_from_meta_lines(meta_lines)
                    validate_versioned_identity(identity)
                except ValueError as exc:
                    raise ValueError(
                        f"Invalid SVCF identity in {str(input_file)!r}: {exc}"
                    ) from exc

                if identity.mode == MODE_MULTI:
                    raise ValueError(
                        f"Input {str(input_file)!r} is a synthesized sample/multi "
                        "SVCF. Merge inputs must be single-caller or per-sample "
                        "caller-merged SVCFs, not previously synthesized "
                        "sample/multi output. Use the original per-sample "
                        "caller-mode SVCF input instead."
                    )

                contig_lengths = _extract_contig_lengths_from_meta_lines(
                    meta_lines, input_file=input_file
                )

                found_chrom_header = True
                header_sample_count = max(0, line.count("\t") - 8)

                try:
                    validate_header_columns_for_identity(
                        identity,
                        header_sample_count,
                    )
                except ValueError as exc:
                    raise ValueError(
                        f"Invalid SVCF header in {str(input_file)!r}: {exc}"
                    ) from exc

                if mode == "sample" and header_sample_count > 1:
                    raise ValueError(
                        f"Input {str(input_file)!r} has {header_sample_count} "
                        "sample columns. `merge --mode sample` currently "
                        "expects one input file per biological sample. "
                        "Split or subset the file to one sample per input "
                        "before merging."
                    )
                continue

            if line.startswith("#") or not line.strip():
                continue

            if not found_chrom_header:
                raise ValueError(
                    f"Input {str(input_file)!r} is missing a #CHROM header "
                    "before its data records."
                )

            fields = line.rstrip("\r\n").split("\t")
            if len(fields) < 10:
                raise ValueError(
                    f"Input {str(input_file)!r} has a malformed data record on "
                    f"line {line_number}: expected at least 10 columns."
                )

            format_field = fields[8]
            try:
                schema_mode = validate_format_for_identity(
                    identity,
                    format_field,
                )
            except ValueError as exc:
                raise ValueError(
                    f"Invalid SVCF schema in {str(input_file)!r} on data line "
                    f"{line_number}: {exc}"
                ) from exc

            if schema_mode is None:
                raise ValueError(
                    f"Input {str(input_file)!r} is not recognized as SVCF on "
                    f"data line {line_number}: FORMAT={format_field!r}. "
                    "Run `octopusv correct` on ordinary caller VCF input first."
                )

            if schema_mode == MODE_MULTI:
                raise ValueError(
                    f"Input {str(input_file)!r} contains the synthesized "
                    f"sample-mode SVCF schema on data line {line_number}. "
                    "Sample-level calls must not be re-used as caller evidence."
                )

            evidence_count = max(0, len(fields) - 9)
            if mode == "caller" and evidence_count > 1:
                raise ValueError(
                    f"Input {str(input_file)!r} contains {evidence_count} "
                    f"evidence/sample columns on data line {line_number}. "
                    "Caller-mode re-merge of multi-evidence SVCF records "
                    "is not supported because nested evidence cannot be "
                    "preserved unambiguously. Use single-evidence SVCF "
                    "inputs for `merge --mode caller`."
                )

    if not found_chrom_header:
        raise ValueError(
            f"Input {str(input_file)!r} is missing the required #CHROM header."
        )

    return contig_lengths


def _preflight_merge_inputs(
    *,
    input_files: list[Path | str],
    labels: list[str],
    mode: str,
    specific: list[Path | str] | None = None,
    expression: str | None = None,
) -> dict[str, str]:
    """Run merge-specific identity, selection, shape, and reference checks."""
    _validate_distinct_input_files(input_files)
    _validate_unique_display_labels(
        labels=labels,
        input_files=input_files,
        mode=mode,
    )
    _validate_display_labels(labels=labels, mode=mode)
    validate_selection_inputs(
        input_files=input_files,
        specific=specific,
        expression=expression,
    )

    merged_contigs: dict[str, str] = {}
    contig_source: dict[str, str] = {}
    for input_file in input_files:
        file_contigs = _preflight_svcf_shape(input_file, mode)
        for contig_id, contig_length in file_contigs.items():
            previous = merged_contigs.get(contig_id)
            if previous is not None and int(previous) != int(contig_length):
                raise ValueError(
                    f"Conflicting ##contig lengths across merge inputs for "
                    f"{contig_id!r}: {previous} in {contig_source[contig_id]!r} "
                    f"versus {contig_length} in {str(input_file)!r}. "
                    "Inputs may use different reference assemblies."
                )
            if previous is None:
                merged_contigs[contig_id] = contig_length
                contig_source[contig_id] = str(input_file)

    return merged_contigs

def _preserve_sample_input_evidence(events) -> None:
    """Preserve exact per-input caller evidence for later sample synthesis.

    The historical function name is kept temporarily because it is internal and
    already used by the merge command/tests, but 1.0 no longer performs a
    first-evidence collapse here. Instead, every parsed sample-mode input record
    receives one runtime-only payload containing:

    * the original caller evidence blocks and explicit source mapping; and
    * the record-level fields needed to build a synthesized sample call after
      core merge/representative selection has finished.

    No public sample field is changed at this stage. In particular, ``SC`` and
    the first evidence block stay untouched so core merge behavior remains
    identical. Runtime-only ``_octopusv_`` keys must never be serialized.
    """
    payload_key = "_octopusv_evidence_payload"

    for event in events:
        raw_columns = list(getattr(event, "raw_sample_columns", []) or [])
        sample_data = getattr(event, "sample", None)
        if not isinstance(sample_data, dict):
            continue

        info = getattr(event, "info", None)
        if not isinstance(info, dict):
            info = {}

        sample_data[payload_key] = {
            "format": getattr(event, "format", ""),
            "blocks": tuple(raw_columns),
            "sources": info.get("SOURCES"),
            "source_ids": info.get("SOURCE_IDS"),
            "record": {
                "id": getattr(event, "sv_id", "."),
                "ref": getattr(event, "ref", "."),
                "alt": getattr(event, "alt", "."),
                "quality": getattr(event, "quality", "."),
                "svtype": info.get("SVTYPE", "."),
                "strand": info.get("STRAND", "."),
                "svlen": info.get("SVLEN", "."),
                "chrom": getattr(event, "chrom", "."),
                "pos": getattr(event, "pos", "."),
                "end_chrom": getattr(event, "end_chrom", info.get("CHR2", ".")),
                "end_pos": getattr(event, "end_pos", info.get("END", ".")),
            },
        }


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
    labels_from_input_list: bool = False,
) -> None:
    """Explain where output caller/sample labels came from."""
    default_labels = _default_labels_from_input_files(input_files)

    if labels_from_input_list:
        _echo(
            "Note: Using caller/sample labels from the TAB-separated second "
            "column of --input-list."
        )
        return

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
    labels_from_input_list: bool,
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
        labels_from_input_list=labels_from_input_list,
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
        with open_text_auto(filename) as f:
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
        input_list: Path | None = typer.Option(
            None,
            "--input-list",
            help=(
                "Text file containing one input path per line, optionally followed "
                "by a TAB and caller/sample label. Blank lines and lines whose first "
                "non-whitespace character is '#' are ignored. Relative paths are "
                "resolved relative to the list file."
            ),
        ),
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

    # Explicit label vectors are positionally bound to inputs. Mixing multiple
    # input mechanisms would make that binding dependent on CLI ordering, so
    # require one mechanism whenever --caller-names/--sample-names is used.
    if caller_names or sample_names:
        input_mechanism_count = sum(
            [bool(input_files), bool(input_option), input_list is not None]
        )
        if input_mechanism_count > 1:
            label_option = "--caller-names" if caller_names else "--sample-names"
            typer.echo(
                f"Error: {label_option} cannot be combined with multiple input "
                "mechanisms. Use only positional inputs, only -i/--input-file, "
                "or only --input-list so labels cannot be bound to the wrong file.",
                err=True,
            )
            raise typer.Exit(code=1)

    # Handle direct inputs while preserving the existing unlabeled CLI ordering rule.
    all_input_files: list[Path] = []
    using_i_option = any(
        arg == "-i" or arg == "--input-file" or arg.startswith("--input-file=")
        for arg in sys.argv[1:]
    )

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

    list_labels: list[str] | None = None
    if input_list is not None:
        try:
            listed_files, list_labels = _read_merge_input_list(input_list)
        except (OSError, ValueError) as exc:
            _echo(f"Error: {exc}")
            raise typer.Exit(code=1)

        if list_labels is not None:
            if all_input_files:
                _echo(
                    "Error: A labeled --input-list cannot be combined with direct "
                    "positional/-i inputs. Put all inputs in the labeled list."
                )
                raise typer.Exit(code=1)
            if caller_names or sample_names:
                _echo(
                    "Error: Labels from --input-list cannot be combined with "
                    "--caller-names or --sample-names."
                )
                raise typer.Exit(code=1)
            all_input_files = listed_files
        else:
            # Unlabeled lists may be combined with direct inputs. Direct inputs
            # remain first; list entries follow in file order.
            all_input_files.extend(listed_files)

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

    # Build name mapper. A labeled input list supplies the complete
    # caller/sample label vector directly.
    name_mapper = None
    try:
        if list_labels is not None:
            name_mapper = NameMapper(
                all_input_files,
                mode=mode,
                custom_names=list_labels,
            )
        elif mode == "caller" and caller_names:
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
        contigs = _preflight_merge_inputs(
            input_files=all_input_files,
            labels=labels,
            mode=mode,
            specific=specific,
            expression=expression,
        )
    except (OSError, ValueError) as e:
        _echo(f"Error: {e}")
        raise typer.Exit(code=1)

    _print_initial_configuration(mode=mode, labels=labels, input_files=all_input_files)

    # Process SV events. Header contigs were already collected during preflight.
    input_filenames = [str(file) for file in all_input_files]

    sv_event_creator = SVCFFileEventCreator(input_filenames)
    sv_event_creator.parse()

    if mode == "sample":
        _preserve_sample_input_evidence(sv_event_creator.events)

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

    # Apply merge strategy. Selection errors are user-facing CLI errors, not
    # internal tracebacks. Most argument-only failures are already caught by
    # preflight before parsing/merging; this guard also covers evaluation-time
    # expression errors.
    try:
        if expression:
            results = sv_merger.get_events_by_expression(expression)
        elif intersect:
            results = sv_merger.get_events_by_source(
                [str(file) for file in all_input_files],
                operation="intersection",
            )
        elif union:
            results = sv_merger.get_events_by_source(
                [str(file) for file in all_input_files],
                operation="union",
            )
        elif specific:
            specific_files = [str(file) for file in specific]
            results = sv_merger.get_events_by_source(
                specific_files,
                operation="specific",
            )
        elif exact_support is not None:
            results = sv_merger.get_events_by_exact_support(exact_support)
        elif min_support is not None or max_support is not None:
            results = sv_merger.get_events_by_support_range(min_support, max_support)
        else:
            raise ValueError(
                "No merge strategy specified. Please use --intersect, --union, --specific, "
                "--min-support, --exact-support, --max-support, or --expression."
            )
    except ValueError as exc:
        _echo(f"Error: {exc}")
        raise typer.Exit(code=1)

    # Write merged results atomically. A writer failure must never leave a
    # truncated file at the requested output path. User-facing schema/I/O
    # failures are reported cleanly rather than as tracebacks.
    try:
        with atomic_output_path(output_file) as temp_output:
            sv_merger.write_results(
                temp_output,
                results,
                contigs,
                mode,
                name_mapper,
                input_filenames,
            )
    except (OSError, ValueError) as exc:
        _echo(f"Error: {exc}")
        raise typer.Exit(code=1)

    _print_success_message(
        output_file=output_file,
        result_count=len(results),
        input_count=len(all_input_files),
        mode=mode,
        labels=labels,
        input_files=all_input_files,
        caller_names=caller_names,
        sample_names=sample_names,
        labels_from_input_list=list_labels is not None,
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
