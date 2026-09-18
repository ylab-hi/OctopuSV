from datetime import datetime

from octopusv.utils.text_io import open_text_auto

from octopusv.utils.svcf_schema import (
    MODE_CALLER,
    MODE_MULTI,
    mode_header,
    version_header,
)


SAFE_GLOBAL_META_KEYS = ("reference", "assembly")
_GLOBAL_META_AUDIT_PREFIX = "OctopuSV_input_"


def _quote_meta_value(value):
    """Quote one value for a structured VCF meta-information field.

    This is serialization only.  The value used for equality decisions remains
    the exact raw header value (apart from line terminators).
    """
    escaped = str(value).replace("\\", "\\\\").replace('"', '\\"')
    return f'"{escaped}"'


def global_meta_audit_lines(meta_lines):
    """Return existing OctopuSV input-reference/assembly audit lines verbatim."""
    prefixes = tuple(
        f"##{_GLOBAL_META_AUDIT_PREFIX}{key}="
        for key in SAFE_GLOBAL_META_KEYS
    )
    result = []
    seen = set()
    for raw_line in meta_lines or []:
        line = str(raw_line).rstrip("\r\n")
        if line.startswith(prefixes) and line not in seen:
            seen.add(line)
            result.append(line)
    return result


def merge_safe_global_meta_lines(meta_sources):
    """Preserve verifiable global metadata without guessing compatibility.

    ``meta_sources`` is an iterable of ``(source_name, meta_lines)`` pairs.
    Only ``##reference`` and ``##assembly`` are considered here.  Missing
    declarations are allowed.  If all observed values for one key are exactly
    identical, that standard header line is preserved.

    If multiple distinct raw strings are observed, OctopuSV cannot determine
    whether those strings describe compatible coordinate systems.  Such a
    difference therefore must not block SV merging and must not be collapsed by
    normalization or first-wins behavior.  Instead, the ambiguous standard
    declaration is omitted, every observed source/value pair is retained in
    OctopuSV audit meta-information, and a warning is emitted.

    Equality is exact.  Values are not case-folded, path-normalized, basename-
    normalized, URI-normalized, or otherwise reinterpreted.  Only ``\r``/``\n``
    line terminators are removed while reading header lines.
    """
    import logging

    observations = {key: [] for key in SAFE_GLOBAL_META_KEYS}

    for source_name, meta_lines in meta_sources:
        source_label = str(source_name) if source_name is not None else "<unknown>"
        for raw_line in meta_lines or []:
            line = str(raw_line).rstrip("\r\n")
            for key in SAFE_GLOBAL_META_KEYS:
                prefix = f"##{key}="
                if not line.startswith(prefix):
                    continue
                value = line[len(prefix):]
                observations[key].append((source_label, value))
                break

    output_lines = []

    for key in SAFE_GLOBAL_META_KEYS:
        items = observations[key]
        if not items:
            continue

        unique_values = {value for _source, value in items}
        if len(unique_values) == 1:
            value = next(iter(unique_values))
            output_lines.append(f"##{key}={value}")
            continue

        # Deduplicate identical source/value pairs but preserve every distinct
        # observation.  Sort only for deterministic serialization; source and
        # value strings themselves are not normalized.
        distinct_items = sorted(set(items), key=lambda item: (item[0], item[1]))
        details = "; ".join(
            f"{source!r}={value!r}"
            for source, value in distinct_items
        )
        logging.warning(
            "Conflicting ##%s metadata strings across inputs; omitting a "
            "single ##%s declaration because OctopuSV cannot verify their "
            "reference compatibility. Preserving all observed values in "
            "##%s%s audit lines. Inputs: %s",
            key,
            key,
            _GLOBAL_META_AUDIT_PREFIX,
            key,
            details,
        )

        audit_key = f"{_GLOBAL_META_AUDIT_PREFIX}{key}"
        for source, value in distinct_items:
            output_lines.append(
                f"##{audit_key}=<Source={_quote_meta_value(source)},"
                f"Value={_quote_meta_value(value)}>"
            )

    return output_lines


def extract_original_header_definitions(input_vcf_file):
    """
    Extract header definitions from original VCF file.
    Returns dictionary with different types of header lines.

    🔴 CHANGED: Also extracts the list of sample names from the #CHROM
    column header. For a single-sample VCF this list has one entry; for a
    multi-sample VCF it has one entry per sample column.
    """
    header_info = {
        'filter_lines': [],
        'info_lines': [],
        'format_lines': [],
        'contig_lines': [],
        'alt_lines': [],
        'other_lines': [],
        # 🔴 CHANGED: default to a single "Sample" column to preserve the
        # legacy behavior when no #CHROM line is encountered.
        'sample_names': ["Sample"],
    }

    is_svaba_output = False

    with open_text_auto(input_vcf_file) as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith("##source=") and "svaba" in line.lower():
                is_svaba_output = True

            # Parse sample names from the #CHROM line, then stop.  Preserve the
            # historical SvABA 13-column compatibility path used by the parser:
            # that layout intentionally selects only the final sample column.
            # Keeping the header selection aligned with the parser prevents a
            # malformed SVCF with four declared samples but one sample block.
            if line.startswith('#CHROM'):
                cols = line.split("\t")
                if len(cols) >= 10:
                    if is_svaba_output and len(cols) == 13:
                        header_info['sample_names'] = [cols[12]]
                    else:
                        header_info['sample_names'] = cols[9:]
                break
            if not line.startswith('##'):
                # Reached a data line without a #CHROM header; stop anyway.
                break

            if line.startswith('##FILTER='):
                header_info['filter_lines'].append(line)
            elif line.startswith('##INFO='):
                header_info['info_lines'].append(line)
            elif line.startswith('##FORMAT='):
                header_info['format_lines'].append(line)
            elif line.startswith('##contig='):
                header_info['contig_lines'].append(line)
            elif line.startswith('##ALT='):
                header_info['alt_lines'].append(line)
            else:
                # Other header lines (fileformat, source, etc.)
                header_info['other_lines'].append(line)

    return header_info


def get_octopus_default_definitions():
    """
    Get OctopuSV default header definitions.
    Returns dictionary with default definitions that OctopuSV adds.
    """
    return {
        'alt_lines': [
            '##ALT=<ID=DEL,Description="Deletion">',
            '##ALT=<ID=INV,Description="Inversion">',
            '##ALT=<ID=INS,Description="Insertion">',
            '##ALT=<ID=DUP,Description="Duplication">',
            '##ALT=<ID=TRA,Description="Translocation">',
            '##ALT=<ID=BND,Description="Breakend">'
        ],
        'info_lines': [
            '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">',
            '##INFO=<ID=CHR2,Number=1,Type=String,Description="Chromosome for end">',
            '##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant described in this record">',
            '##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Difference in length between REF and ALT alleles">',
            '##INFO=<ID=SUPPORT,Number=1,Type=Integer,Description="Number of pieces of evidence supporting the variant">',
            '##INFO=<ID=SVMETHOD,Number=1,Type=String,Description="The software used to identify the SV">',
            '##INFO=<ID=RTID,Number=1,Type=String,Description="Associated ID for reciprocal translocations if available">',
            '##INFO=<ID=AF,Number=1,Type=Float,Description="Allele Frequency">',
            '##INFO=<ID=STRAND,Number=1,Type=String,Description="Strand orientation of the SV">',
            '##INFO=<ID=RNAMES,Number=.,Type=String,Description="Supporting read names">' ,
            '##INFO=<ID=SOURCES,Number=.,Type=String,Description="Source caller/sample labels supporting this merged SV record">',
            '##INFO=<ID=SOURCE_IDS,Number=.,Type=String,Description="Original IDs of merged SVs from different callers or samples">'
        ],
        'format_lines': [
            '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
            '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">',
            '##FORMAT=<ID=LN,Number=1,Type=Integer,Description="Length of SV">',
            '##FORMAT=<ID=ST,Number=1,Type=String,Description="Strand orientation of SV (e.g., +, -, -+, ++)">',
            '##FORMAT=<ID=QV,Number=1,Type=Float,Description="Quality value">',
            '##FORMAT=<ID=TY,Number=1,Type=String,Description="Type of SV (e.g., TRA, DEL, INS)">',
            '##FORMAT=<ID=ID,Number=1,Type=String,Description="Unique identifier for the SV">',
            '##FORMAT=<ID=SC,Number=1,Type=String,Description="Source from which SV was identified">',
            '##FORMAT=<ID=REF,Number=1,Type=String,Description="Reference allele sequence">',
            '##FORMAT=<ID=ALT,Number=1,Type=String,Description="Alternate allele sequence">',
            '##FORMAT=<ID=CO,Number=1,Type=String,Description="Coordinate information of the SV">'
        ],
        'filter_lines': [
            '##FILTER=<ID=PASS,Description="All filters passed, variant is most likely true">'
        ]
    }


def extract_id_from_header_line(line, prefix):
    """
    Extract ID from header line.
    Example: ##INFO=<ID=SVTYPE,... -> returns "SVTYPE"
    """
    try:
        if f'##{prefix}=<ID=' in line:
            return line.split('ID=')[1].split(',')[0]
    except IndexError:
        pass
    return None


def merge_header_definitions(original_headers, octopus_defaults):
    """Merge input headers while keeping SVCF canonical definitions authoritative.

    Preserve the historical input-header ordering as much as possible. When an
    input definition reuses an ID owned by OctopuSV, emit the canonical SVCF
    definition at that same position instead of the conflicting input line.
    Canonical definitions absent from the input are appended afterward in their
    normal OctopuSV order. This avoids broad header reordering while preventing
    raw/legacy definitions from overriding the SVCF contract.
    """
    merged = {
        "filter_lines": [],
        "info_lines": [],
        "format_lines": [],
        "alt_lines": [],
    }

    categories = (
        ("filter_lines", "FILTER"),
        ("info_lines", "INFO"),
        ("format_lines", "FORMAT"),
        ("alt_lines", "ALT"),
    )

    for key, prefix in categories:
        default_lines = list(octopus_defaults.get(key, []))
        default_by_id = {}
        for line in default_lines:
            item_id = extract_id_from_header_line(line, prefix)
            if item_id is not None:
                default_by_id[item_id] = line

        emitted_ids = set()
        emitted_raw_lines = set()

        # Walk input definitions first so non-conflicting caller metadata keeps
        # its historical relative ordering. A conflicting SVCF-owned ID is
        # replaced in place by the canonical definition.
        for line in original_headers.get(key, []):
            item_id = extract_id_from_header_line(line, prefix)
            if item_id is not None:
                if item_id in emitted_ids:
                    continue
                if item_id in default_by_id:
                    merged[key].append(default_by_id[item_id])
                else:
                    merged[key].append(line)
                emitted_ids.add(item_id)
                continue

            if line in emitted_raw_lines:
                continue
            merged[key].append(line)
            emitted_raw_lines.add(line)

        # Match historical behavior for missing OctopuSV definitions: append
        # them after the preserved input definitions, in canonical order.
        for line in default_lines:
            item_id = extract_id_from_header_line(line, prefix)
            if item_id is not None and item_id in emitted_ids:
                continue
            if item_id is None and line in emitted_raw_lines:
                continue
            merged[key].append(line)
            if item_id is not None:
                emitted_ids.add(item_id)
            else:
                emitted_raw_lines.add(line)

    return merged


def generate_sv_header(contig_lines, input_vcf_file=None, extra_meta_lines=None):
    """
    Generate SVCF file header lines according to SVCF specification.
    If input_vcf_file is provided, extract and preserve original header definitions.

    SVCF identity follows the data model, not merely the command that wrote
    the file:

    * a single-sample ``correct`` output is caller-evidence SVCF 1.1 and
      explicitly declares ``SVCFVersion=1.1`` plus ``OctopuSV_mode=caller``;
    * a multi-sample raw-VCF ``correct`` output keeps the historical
      unversioned ``OctopuSV_mode=multi`` marker.  Its columns are biological
      samples carrying caller-style evidence blocks, so it is neither the
      caller-evidence matrix nor the synthesized multi-sample model defined by
      SVCF 1.1.

    The #CHROM line always preserves the original input sample names (falling
    back to ``Sample`` only when no input file is supplied).
    """
    current_time_str = datetime.now().strftime("%Y-%m-%d|%I:%M:%S%p|%Z")

    # Basic header lines
    basic_header = [
        "##fileformat=VCFv4.2",
        f"##fileDate={current_time_str}",
        "##source=OctopuSV",
        "##OctopuSV_WARNING=This is SVCF format. Use 'octopusv svcf2vcf' to change back to standard VCF format before bcftools/vcftools"
    ]

    # Get OctopuSV default definitions
    octopus_defaults = get_octopus_default_definitions()

    # 🔴 CHANGED: defaults that get overridden when we have an input file.
    sample_names = ["Sample"]
    is_multi_sample = False

    # Preserve only narrowly approved global metadata. Generic ``other_lines``
    # may contain stale caller/tool identity and are deliberately not copied.
    safe_global_meta_lines = []

    # If input VCF file is provided, extract original definitions
    if input_vcf_file:
        original_headers = extract_original_header_definitions(input_vcf_file)
        safe_global_meta_lines = merge_safe_global_meta_lines(
            [(str(input_vcf_file), original_headers.get('other_lines', []))]
        )
        # Use original contig lines if available, otherwise use provided ones
        if original_headers['contig_lines']:
            contig_lines = original_headers['contig_lines']
        # 🔴 CHANGED: take sample names from the original #CHROM line.
        sample_names = original_headers.get('sample_names', ["Sample"])
        is_multi_sample = len(sample_names) > 1
        merged_definitions = merge_header_definitions(original_headers, octopus_defaults)
    else:
        # Use only OctopuSV defaults
        merged_definitions = octopus_defaults

    # Identity is explicit only when this output genuinely satisfies one of
    # the SVCF 1.1 data models.  Single-sample correct output is one caller
    # observation/evidence column and therefore qualifies as caller mode.
    # Historical multi-sample correct output is intentionally left
    # unversioned: its columns are biological samples, but its blocks still use
    # caller evidence FORMAT rather than synthesized UC/UV sample calls.
    final_header = list(basic_header)
    final_header.extend(safe_global_meta_lines)
    if extra_meta_lines:
        final_header.extend(extra_meta_lines)
    if is_multi_sample:
        final_header.append(mode_header(MODE_MULTI))
    else:
        final_header.append(version_header())
        final_header.append(mode_header(MODE_CALLER))

    final_header.extend(contig_lines)

    # Add definitions in standard order
    final_header.extend(merged_definitions['alt_lines'])
    final_header.extend(merged_definitions['info_lines'])
    final_header.extend(merged_definitions['filter_lines'])
    final_header.extend(merged_definitions['format_lines'])

    # 🔴 CHANGED: build the #CHROM line dynamically so it contains the
    # correct number of sample columns. Single-sample VCFs end up with
    # exactly one trailing column (matching the legacy behavior, but now
    # carrying the original sample name instead of the placeholder "Sample").
    chrom_cols = ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"]
    chrom_cols.extend(sample_names)
    final_header.append("\t".join(chrom_cols))

    return final_header


def write_sv_vcf(contig_lines, events, output_file, input_vcf_file=None, extra_meta_lines=None):
    """
    Write SV events to VCF file with proper header definitions.
    If input_vcf_file is provided, preserve original header definitions.

    Records are written in stable genomic order using the exact ``##contig``
    declaration order followed by POS. Sorting is output-only and never
    mutates event contents.
    """
    from octopusv.utils.svcf_sort import (
        contig_order_from_meta_lines,
        sort_events_for_output,
    )

    sv_header = generate_sv_header(contig_lines, input_vcf_file, extra_meta_lines=extra_meta_lines)
    ordered_events = sort_events_for_output(
        events,
        contig_order_from_meta_lines(contig_lines),
    )
    with open(output_file, "w") as f:
        for header in sv_header:
            f.write(header + "\n")
        for event in ordered_events:
            f.write(str(event) + "\n")
