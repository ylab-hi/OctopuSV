import logging
import re
import sys

from octopusv.converter.base import get_bnd_pattern
from octopusv.sv import SVEvent
from octopusv.utils.text_io import open_text_auto
from octopusv.utils.svcf_schema import validate_source_label


_SINGLE_BREAKEND_RE = re.compile(r"^(?:[ACGTNacgtn]+\.|\.[ACGTNacgtn]+)$")
_SPAN_SVTYPES = {"DEL", "DUP", "INV"}
_NON_BND_SVTYPES = {"DEL", "DUP", "INV", "INS", "TRA"}
_LEGAL_SVTYPES = _NON_BND_SVTYPES | {"BND"}
_REFERENCE_SYMBOLIC_ALTS = {"<NON_REF>", "<*>"}
_MISSING_EXPLICIT_GT = {"", ".", "./.", ".|."}


def _synthesize_site_only_sample(format_field: str | None) -> tuple[str, str]:
    """Return a FORMAT/sample pair for a site-only VCF record.

    A record that exists in a caller VCF establishes carrier presence, but a
    missing sample column provides no zygosity information.  Represent exactly
    that knowledge as ``GT=1/.`` rather than inventing heterozygosity.

    Nine-column VCFs may already declare FORMAT keys despite having no sample
    column.  Preserve those keys, add GT when absent, and create placeholders
    aligned to the resulting FORMAT so downstream parsing cannot shift fields.
    """
    raw_format = str(format_field or "").strip()
    if raw_format in {"", "."}:
        return "GT", "1/."

    format_keys = raw_format.split(":")
    if "GT" not in format_keys:
        format_keys = ["GT", *format_keys]

    sample_values = []
    for key in format_keys:
        if key == "GT":
            sample_values.append("1/.")
        elif key == "AD":
            # AD is Number=R in the SVCF output; an unknown diploid ref/alt
            # depth therefore remains the explicit two-value missing form.
            sample_values.append(".,.")
        else:
            sample_values.append(".")

    return ":".join(format_keys), ":".join(sample_values)


def _explicit_gt_value(format_field: str, sample_value: str) -> str | None:
    """Return an explicitly supplied GT, or ``None`` when FORMAT has no GT.

    This helper is used only for the user-facing all-missing warning.  Site-only
    records synthesized by OctopuSV never call it, so ``GT=1/.`` created above
    is not mistaken for caller-supplied genotype evidence.
    """
    format_keys = str(format_field or "").split(":")
    if "GT" not in format_keys:
        return None

    gt_index = format_keys.index("GT")
    sample_parts = str(sample_value or "").split(":")
    if gt_index >= len(sample_parts):
        return "."
    return sample_parts[gt_index].strip()


def _is_single_breakend_alt(alt: str) -> bool:
    """Return True for VCF one-ended breakends such as ``N.`` or ``.N``."""
    return bool(_SINGLE_BREAKEND_RE.fullmatch(alt))


def _looks_structural_without_svtype(event) -> bool:
    """Return True when a record clearly looks like an SV despite missing SVTYPE.

    OctopuSV historically skipped records without ``SVTYPE`` so mixed gVCF/
    reference-block inputs remained usable.  Keep that compatibility for
    ordinary non-SV alleles and the common reference placeholders, but do not
    silently discard records that are unmistakably structural variants.
    """
    alt = str(event.alt or "")
    alt_upper = alt.upper()

    # SVLEN and CHR2 are strong structural-variant signals.  END alone is not:
    # gVCF reference blocks routinely carry END without representing an SV.
    if "SVLEN" in event.info or "CHR2" in event.info:
        return True

    if alt_upper in _REFERENCE_SYMBOLIC_ALTS or alt in {".", "*"}:
        return False

    if _is_single_breakend_alt(alt) or "[" in alt or "]" in alt:
        return True

    # Any other symbolic allele is structural-like in the context of
    # ``octopusv correct``.  We deliberately do not guess which SVTYPE it is.
    if alt.startswith("<") and alt.endswith(">"):
        return True

    return False


def _canonicalize_and_validate_svtype(event, line_number: int) -> str:
    """Return the canonical SVTYPE that ``correct`` is allowed to emit.

    Preserve the historical non-BND ``TYPE:subtype`` simplification used by
    ``NonBNDConverter`` while refusing unsupported final types instead of
    writing an SVCF that the validator will later reject.
    """
    raw_svtype = str(event.info.get("SVTYPE", ""))
    simplified = raw_svtype.split(":", 1)[0]

    if raw_svtype == "BND":
        return "BND"

    if simplified in _NON_BND_SVTYPES:
        event.info["SVTYPE"] = simplified
        return simplified

    raise ValueError(
        f"Unsupported SVTYPE at line {line_number} "
        f"({event.id} {event.chrom}:{event.pos}): SVTYPE={raw_svtype!r}. "
        f"OctopuSV correct can emit only {sorted(_LEGAL_SVTYPES)}."
    )


def _normalize_span_end(event, line_number: int) -> None:
    """Ensure DEL/DUP/INV have a trustworthy numeric END.

    A missing END can be reconstructed exactly when an integer SVLEN is
    available.  An explicitly malformed or contradictory END is never
    overwritten: that would replace caller-provided information with a guess.
    """
    if event.info.get("SVTYPE") not in _SPAN_SVTYPES:
        return

    end = event.info.get("END")
    if end in (None, "", "."):
        svlen = event.info.get("SVLEN")
        try:
            svlen_int = int(str(svlen))
        except (TypeError, ValueError):
            raise ValueError(
                f"{event.info['SVTYPE']} record at line {line_number} "
                f"({event.id} {event.chrom}:{event.pos}) has no numeric END and "
                f"SVLEN={svlen!r} cannot be used to derive one."
            ) from None

        event.info["END"] = str(event.pos + abs(svlen_int))
        return

    try:
        end_int = int(str(end))
    except (TypeError, ValueError):
        raise ValueError(
            f"Malformed END at line {line_number} "
            f"({event.id} {event.chrom}:{event.pos}): END={end!r}. "
            "OctopuSV will not replace an explicitly malformed END with a guessed coordinate."
        ) from None

    if end_int < event.pos:
        raise ValueError(
            f"Invalid END at line {line_number} "
            f"({event.id} {event.chrom}:{event.pos}): END={end_int} < POS={event.pos}."
        )


def _normalize_and_validate_tra(event, line_number: int) -> None:
    """Preserve TRA with two known breakpoints without inventing orientation.

    Symbolic ``<TRA>`` requires explicit CHR2/END.  Bracket-form TRA may safely
    supply missing CHR2/END from ALT because those coordinates are encoded in
    the record itself; explicit conflicts are rejected rather than guessed.
    """
    if event.info.get("SVTYPE") != "TRA":
        return

    chr2 = event.info.get("CHR2")
    end = event.info.get("END")

    if event.alt == "<TRA>":
        if chr2 in (None, "", "."):
            raise ValueError(
                f"TRA record at line {line_number} "
                f"({event.id} {event.chrom}:{event.pos}) uses ALT=<TRA> but has no CHR2. "
                "Two breakpoint coordinates are required; orientation may remain unknown."
            )
        if end in (None, "", "."):
            raise ValueError(
                f"TRA record at line {line_number} "
                f"({event.id} {event.chrom}:{event.pos}) uses ALT=<TRA> but has no numeric END. "
                "Two breakpoint coordinates are required; orientation may remain unknown."
            )
        try:
            int(str(end))
        except (TypeError, ValueError):
            raise ValueError(
                f"TRA record at line {line_number} "
                f"({event.id} {event.chrom}:{event.pos}) has non-numeric END={end!r}."
            ) from None
        return

    chrom_alt, pos_alt = _paired_bnd_remote_target(event.alt)
    if chrom_alt is None or pos_alt is None:
        raise ValueError(
            f"Malformed or unsupported TRA ALT at line {line_number} "
            f"({event.id} {event.chrom}:{event.pos} ALT={event.alt}). "
            "TRA must use either symbolic <TRA> with CHR2/END or a valid bracket ALT."
        )

    if chr2 in (None, "", "."):
        event.info["CHR2"] = chrom_alt
    elif str(chr2) != chrom_alt:
        raise ValueError(
            f"Conflicting TRA mate chromosome at line {line_number} "
            f"({event.id} {event.chrom}:{event.pos}): ALT points to {chrom_alt}, "
            f"but CHR2={chr2}."
        )

    if end in (None, "", "."):
        event.info["END"] = str(pos_alt)
    else:
        try:
            end_int = int(str(end))
        except (TypeError, ValueError):
            raise ValueError(
                f"TRA record at line {line_number} "
                f"({event.id} {event.chrom}:{event.pos}) has non-numeric END={end!r}."
            ) from None
        if end_int != pos_alt:
            raise ValueError(
                f"Conflicting TRA mate position at line {line_number} "
                f"({event.id} {event.chrom}:{event.pos}): ALT points to {pos_alt}, "
                f"but END={end_int}."
            )


def _paired_bnd_remote_target(alt: str):
    """Extract ``(chrom, pos)`` from a paired BND ALT for validation.

    This helper deliberately does not broaden the accepted BND grammar.  It
    parses the final ``:pos`` delimiter from the right only so OctopuSV can
    distinguish a currently unsupported colon-containing remote contig from a
    generic malformed ALT and fail explicitly instead of silently dropping the
    record.
    """
    match = re.search(r"[\[\]]([^\[\]]+)[\[\]]", alt)
    if not match:
        return None, None

    remote = match.group(1)
    if ":" not in remote:
        return None, None

    chrom, pos_text = remote.rsplit(":", 1)
    if not chrom or not pos_text.isdigit():
        return None, None
    return chrom, int(pos_text)


def _validate_bnd_alt(event, line_number: int) -> str:
    """Classify a raw BND ALT before it enters the correction pipeline.

    ``paired`` means one of OctopuSV's existing four paired-breakend forms;
    ``single`` means a true one-ended VCF breakend.  Other forms are rejected
    explicitly so malformed input or parser limitations cannot become silent
    record loss.
    """
    alt = event.alt

    if _is_single_breakend_alt(alt):
        return "single"

    pattern = get_bnd_pattern(alt)
    if pattern is None:
        raise ValueError(
            f"Malformed or unsupported BND ALT at line {line_number} "
            f"({event.id} {event.chrom}:{event.pos} ALT={alt}). "
            "OctopuSV correct will not silently discard an unparseable breakend."
        )

    chrom_alt, pos_alt = _paired_bnd_remote_target(alt)
    if chrom_alt is None or pos_alt is None:
        raise ValueError(
            f"Malformed or unsupported BND ALT at line {line_number} "
            f"({event.id} {event.chrom}:{event.pos} ALT={alt})."
        )

    if ":" in chrom_alt:
        raise ValueError(
            f"Unsupported BND remote contig at line {line_number} "
            f"({event.id} {event.chrom}:{event.pos} ALT={alt}): "
            "contig names containing ':' are not yet supported by the full "
            "OctopuSV SVCF/BND pipeline. Filter those contigs before running "
            "octopusv correct, or use a reference naming scheme without ':'."
        )

    return "paired"


def _validate_svcf_contig_names(event, line_number: int) -> None:
    """Reject contig names that the current SVCF evidence encoding cannot round-trip.

    SVCF caller evidence blocks use ``:`` as a field delimiter and the current
    shared parser has no escaping layer for ``CO``.  A local ``CHROM`` or an
    explicit ``INFO/CHR2`` containing ``:`` would therefore serialize
    successfully but be silently misparsed downstream.  Until full-stack
    support is added, fail at the raw-VCF boundary instead of writing a
    corrupted SVCF representation.
    """
    chr2 = event.info.get("CHR2")
    offending = None
    label = None

    if ":" in event.chrom:
        offending = event.chrom
        label = "CHROM"
    elif chr2 not in (None, ".") and ":" in str(chr2):
        offending = str(chr2)
        label = "CHR2"

    if offending is not None:
        raise ValueError(
            f"Unsupported contig name at line {line_number} "
            f"({event.id} {event.chrom}:{event.pos}): {label}={offending}. "
            "The current OctopuSV SVCF pipeline does not support contig names containing ':' "
            "because they cannot currently be serialized "
            "and parsed losslessly. Filter those contigs before running "
            "octopusv correct, or use a reference naming scheme without ':'."
        )


def is_same_chr_bnd(event):
    """Check if the POS and ALT of an event are on the same chromosome."""
    if event.is_BND():
        split_result = re.split(r"[\[\]:]", event.alt)
        if len(split_result) != 4:
            return False
        chrom_alt, _ = split_result[1:3]
        return event.chrom == chrom_alt

    return False  # For non-BND, we won't categorize them as same_chr_bnd or diff_chr_bnd events


def check_vcf_format(vcf_file_path):
    """Check the format of a VCF file in a streaming pass.

    Plain-text and gzip/bgzip-compressed VCF inputs are accepted.  A real
    ``#CHROM`` column header is required, and every data row must have the
    same number of tab-delimited columns declared by that header.  This keeps
    simplified 8/9-column inputs and standard single-/multi-sample VCFs
    supported while rejecting truncated or structurally inconsistent rows.
    """
    chrom_header_fields = None

    with open_text_auto(vcf_file_path) as f:
        for line_number, line in enumerate(f, start=1):
            if line.startswith("##"):
                continue

            if line.startswith("#CHROM"):
                if chrom_header_fields is not None:
                    logging.error(
                        f"Invalid VCF format at line {line_number}: multiple #CHROM header lines found.",
                    )
                    sys.exit(1)

                chrom_header_fields = line.rstrip("\r\n").split("\t")
                if len(chrom_header_fields) < 8:
                    logging.error(
                        "Invalid VCF format. #CHROM header must contain at least the 8 fixed VCF columns.",
                    )
                    sys.exit(1)
                continue

            if line.startswith("#"):
                # Preserve tolerance for non-standard comment/header lines,
                # but they do not satisfy the required VCF #CHROM header.
                continue

            if not line.strip():
                continue

            if chrom_header_fields is None:
                logging.error(
                    f"Invalid VCF format at line {line_number}: data record encountered before a #CHROM header.",
                )
                sys.exit(1)

            if " " in line:
                logging.error(
                    f"Invalid VCF format at line {line_number}. Non-header lines should not contain spaces.",
                )
                sys.exit(1)

            fields = line.rstrip("\r\n").split("\t")

            if len(fields) < 8:
                logging.error(
                    f"Invalid VCF format at line {line_number}. Expected at least 8 fields, but got {len(fields)}",
                )
                sys.exit(1)

            if len(fields) != len(chrom_header_fields):
                logging.error(
                    f"Invalid VCF format at line {line_number}: #CHROM declares "
                    f"{len(chrom_header_fields)} columns, but the record has {len(fields)}.",
                )
                sys.exit(1)

            try:
                int(fields[1])
            except ValueError:
                logging.error(
                    f"Invalid VCF format at line {line_number}. Position (field 2) should be a number, but got {fields[1]}",
                )
                sys.exit(1)

            if fields[5] != ".":
                try:
                    float(fields[5])
                except ValueError:
                    logging.error(
                        f"Invalid VCF format at line {line_number}. Quality score (field 6) should be a number or '.', but got {fields[5]}",
                    )
                    sys.exit(1)

    if chrom_header_fields is None:
        logging.error(
            "Invalid VCF format. The file must contain a #CHROM column header line.",
        )
        sys.exit(1)


def parse_vcf(vcf_file_path, *, skip_single_breakends=False, parse_stats=None):
    """Parse VCF into SVEvent lists without silently losing BND records.

    The return shape remains backward-compatible: ``(contigs, same_chr_bnd,
    diff_chr_bnd, non_bnd)``.  ``parse_stats`` is an optional mutable mapping
    used by ``octopusv correct`` to report an explicitly requested lossy skip
    of true single-breakends without changing this long-standing API.
    """
    check_vcf_format(vcf_file_path)
    same_chr_bnd_events = []
    diff_chr_bnd_events = []
    non_bnd_events = []
    contig_lines = []
    is_svaba_output = False
    source_info = "."
    single_breakend_count = 0
    single_breakend_examples = []
    skipped_non_sv_count = 0
    skipped_non_sv_examples = []
    recognized_sv_count = 0
    output_sample_count = None
    explicit_gt_count = 0
    explicit_missing_gt_count = 0
    data_record_count = 0

    with open_text_auto(vcf_file_path) as f:
        for line_number, line in enumerate(f, start=1):
            if line.startswith("##source="):
                source_info = line.split("=", 1)[1].split(" ")[0].strip()
                if source_info not in {"", "."}:
                    try:
                        validate_source_label(source_info)
                    except ValueError as exc:
                        raise ValueError(
                            f"Invalid VCF source label {source_info!r}: {exc}"
                        ) from exc
                if "svaba" in line.lower():
                    is_svaba_output = True
            elif line.startswith("##contig"):
                contig_lines.append(line.strip())
            elif line.startswith("#CHROM"):
                header_fields = line.rstrip("\r\n").split("\t")
                if is_svaba_output and len(header_fields) == 13:
                    output_sample_count = 1
                elif len(header_fields) >= 10:
                    output_sample_count = len(header_fields) - 9
                else:
                    # Site-only 8/9-column VCFs are represented by one
                    # synthesized caller evidence column in correct output.
                    output_sample_count = 1
            elif not line.startswith("#"):
                if not line.strip():
                    continue
                data_record_count += 1
                fields = line.strip().split("\t")

                has_explicit_sample_columns = False
                if len(fields) == 8:
                    synthesized_format, synthesized_sample = _synthesize_site_only_sample(None)
                    core_fields = fields[:8] + [synthesized_format]
                    sample_fields = [synthesized_sample]
                elif is_svaba_output and len(fields) == 13:
                    # Preserve the historical SvABA path exactly: column 8 is
                    # FORMAT and column 12 is the one selected sample.
                    core_fields = fields[:8] + [fields[8]]
                    sample_fields = [fields[12]]
                    has_explicit_sample_columns = True
                elif len(fields) >= 10:
                    core_fields = fields[:9]
                    sample_fields = fields[9:]
                    has_explicit_sample_columns = True
                elif len(fields) == 9:
                    synthesized_format, synthesized_sample = _synthesize_site_only_sample(fields[8])
                    core_fields = fields[:8] + [synthesized_format]
                    sample_fields = [synthesized_sample]
                else:
                    continue

                event = SVEvent(
                    *core_fields,
                    sample=sample_fields[0],
                    samples=sample_fields,
                )
                event.source = source_info

                # Preserve compatibility with mixed gVCF/reference-block and
                # ordinary non-SV records, but never silently discard a record
                # that clearly encodes a structural variant while omitting the
                # required SVTYPE label.
                if event.info.get("SVTYPE") in (None, "", "."):
                    if _looks_structural_without_svtype(event):
                        raise ValueError(
                            f"Structural-variant-like record at line {line_number} "
                            f"({event.id} {event.chrom}:{event.pos} ALT={event.alt}) "
                            "is missing SVTYPE. OctopuSV will not guess the event type."
                        )
                    skipped_non_sv_count += 1
                    if len(skipped_non_sv_examples) < 5:
                        skipped_non_sv_examples.append(
                            f"{event.id} {event.chrom}:{event.pos} ALT={event.alt}"
                        )
                    continue

                if event.info.get("SVTYPE") == "CNV":
                    if "<DEL>" in event.alt.upper():
                        event.info["SVTYPE"] = "DEL"
                    elif "<DUP>" in event.alt.upper():
                        event.info["SVTYPE"] = "DUP"
                    elif "LOSS" in event.id.upper():
                        event.info["SVTYPE"] = "DEL"
                    elif "GAIN" in event.id.upper():
                        event.info["SVTYPE"] = "DUP"
                    elif "SVLEN" in event.info:
                        try:
                            svlen = int(event.info["SVLEN"])
                            if svlen < 0:
                                event.info["SVTYPE"] = "DEL"
                            elif svlen > 0:
                                event.info["SVTYPE"] = "DUP"
                        except ValueError:
                            pass
                    if event.info.get("SVTYPE") == "CNV":
                        raise ValueError(
                            f"Could not resolve CNV record at line {line_number} "
                            f"({event.id} {event.chrom}:{event.pos}) to DEL or DUP. "
                            "OctopuSV will not emit unsupported SVTYPE=CNV."
                        )

                _canonicalize_and_validate_svtype(event, line_number)
                _normalize_span_end(event, line_number)
                _normalize_and_validate_tra(event, line_number)

                _validate_svcf_contig_names(event, line_number)

                recognized_sv_count += 1
                if output_sample_count is None:
                    output_sample_count = len(event.samples)

                if has_explicit_sample_columns:
                    for raw_sample in event.samples:
                        gt_value = _explicit_gt_value(event.format, raw_sample)
                        if gt_value is None:
                            continue
                        explicit_gt_count += 1
                        if gt_value in _MISSING_EXPLICIT_GT:
                            explicit_missing_gt_count += 1

                if event.is_BND():
                    bnd_kind = _validate_bnd_alt(event, line_number)
                    if bnd_kind == "single":
                        single_breakend_count += 1
                        if len(single_breakend_examples) < 5:
                            single_breakend_examples.append(
                                f"{event.id} {event.chrom}:{event.pos} ALT={event.alt}"
                            )
                        # A true one-ended breakend has no remote coordinate,
                        # so it cannot satisfy the locked SVCF 1.1 BND contract.
                        continue

                    if is_same_chr_bnd(event):
                        same_chr_bnd_events.append(event)
                    else:
                        diff_chr_bnd_events.append(event)
                else:
                    non_bnd_events.append(event)

    if parse_stats is not None:
        parse_stats["single_breakends"] = single_breakend_count
        parse_stats["single_breakend_examples"] = list(single_breakend_examples)
        parse_stats["skipped_non_sv_records"] = skipped_non_sv_count
        parse_stats["skipped_non_sv_examples"] = list(skipped_non_sv_examples)
        parse_stats["recognized_sv_records"] = recognized_sv_count
        parse_stats["output_sample_count"] = output_sample_count
        parse_stats["explicit_gt_count"] = explicit_gt_count
        parse_stats["explicit_missing_gt_count"] = explicit_missing_gt_count
        parse_stats["data_records"] = data_record_count

    if single_breakend_count and not skip_single_breakends:
        examples = "; ".join(single_breakend_examples)
        raise ValueError(
            f"Found {single_breakend_count} true single-breakend BND record(s) "
            "without a remote breakpoint coordinate. These records cannot be "
            "represented losslessly in SVCF 1.1. "
            f"Examples: {examples}. "
            "Re-run with --skip-single-breakends only if you explicitly want "
            "OctopuSV to omit these one-ended breakends."
        )

    if recognized_sv_count == 0 and data_record_count > 0:
        examples = "; ".join(skipped_non_sv_examples)
        detail = f" Examples: {examples}." if examples else ""
        raise ValueError(
            "Input contained data records but no structural-variant records that "
            "OctopuSV can correct. "
            f"Skipped {skipped_non_sv_count} non-SV/reference record(s) without SVTYPE."
            f"{detail}"
        )

    return contig_lines, same_chr_bnd_events, diff_chr_bnd_events, non_bnd_events
