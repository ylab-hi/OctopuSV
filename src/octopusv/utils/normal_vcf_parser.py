import logging
import re
import sys

from octopusv.converter.base import get_bnd_pattern
from octopusv.sv import SVEvent
from octopusv.utils.text_io import open_text_auto


_SINGLE_BREAKEND_RE = re.compile(r"^(?:[ACGTNacgtn]+\.|\.[ACGTNacgtn]+)$")


def _is_single_breakend_alt(alt: str) -> bool:
    """Return True for VCF one-ended breakends such as ``N.`` or ``.N``."""
    return bool(_SINGLE_BREAKEND_RE.fullmatch(alt))


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

    with open_text_auto(vcf_file_path) as f:
        for line_number, line in enumerate(f, start=1):
            if line.startswith("##source="):
                source_info = line.split("=")[1].split(" ")[0].strip()
                if "svaba" in line.lower():
                    is_svaba_output = True
            elif line.startswith("##contig"):
                contig_lines.append(line.strip())
            elif not line.startswith("#"):
                fields = line.strip().split("\t")

                if len(fields) == 8:
                    core_fields = fields[:8] + ["GT"]
                    sample_fields = ["0/1"]
                elif is_svaba_output and len(fields) == 13:
                    # Preserve the historical SvABA path exactly: column 8 is
                    # FORMAT and column 12 is the one selected sample.
                    core_fields = fields[:8] + [fields[8]]
                    sample_fields = [fields[12]]
                elif len(fields) >= 10:
                    core_fields = fields[:9]
                    sample_fields = fields[9:]
                elif len(fields) == 9:
                    core_fields = fields[:9]
                    sample_fields = ["0/1"]
                else:
                    continue

                event = SVEvent(
                    *core_fields,
                    sample=sample_fields[0],
                    samples=sample_fields,
                )
                event.source = source_info

                # Skip non-variant records (e.g. Dragen REF regions without SVTYPE).
                if "SVTYPE" not in event.info:
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
                        logging.warning(f"Could not determine specific type for CNV event {event.id}, skipping")
                        continue

                _validate_svcf_contig_names(event, line_number)

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

    return contig_lines, same_chr_bnd_events, diff_chr_bnd_events, non_bnd_events
