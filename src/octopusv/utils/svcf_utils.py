from datetime import datetime


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

    with open(input_vcf_file, 'r') as f:
        for line in f:
            line = line.rstrip("\n")
            # 🔴 CHANGED: parse sample names from the #CHROM line, then stop.
            if line.startswith('#CHROM'):
                cols = line.split("\t")
                if len(cols) >= 10:
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
            '##FORMAT=<ID=QV,Number=1,Type=Integer,Description="Quality value">',
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
    """
    Merge original VCF header definitions with OctopuSV defaults.
    Original definitions take priority to avoid conflicts.
    """
    merged = {
        'filter_lines': [],
        'info_lines': [],
        'format_lines': [],
        'alt_lines': []
    }

    # Track existing IDs to avoid duplicates
    existing_filter_ids = set()
    existing_info_ids = set()
    existing_format_ids = set()
    existing_alt_ids = set()

    # Add original FILTER definitions
    for line in original_headers.get('filter_lines', []):
        filter_id = extract_id_from_header_line(line, 'FILTER')
        if filter_id:
            existing_filter_ids.add(filter_id)
        merged['filter_lines'].append(line)

    # Add original INFO definitions
    for line in original_headers.get('info_lines', []):
        info_id = extract_id_from_header_line(line, 'INFO')
        if info_id:
            existing_info_ids.add(info_id)
        merged['info_lines'].append(line)

    # Add original FORMAT definitions
    for line in original_headers.get('format_lines', []):
        format_id = extract_id_from_header_line(line, 'FORMAT')
        if format_id:
            existing_format_ids.add(format_id)
        merged['format_lines'].append(line)

    # Add original ALT definitions
    for line in original_headers.get('alt_lines', []):
        alt_id = extract_id_from_header_line(line, 'ALT')
        if alt_id:
            existing_alt_ids.add(alt_id)
        merged['alt_lines'].append(line)

    # Add OctopuSV defaults that don't conflict with originals
    for line in octopus_defaults.get('filter_lines', []):
        filter_id = extract_id_from_header_line(line, 'FILTER')
        if filter_id and filter_id not in existing_filter_ids:
            merged['filter_lines'].append(line)

    for line in octopus_defaults.get('info_lines', []):
        info_id = extract_id_from_header_line(line, 'INFO')
        if info_id and info_id not in existing_info_ids:
            merged['info_lines'].append(line)

    for line in octopus_defaults.get('format_lines', []):
        format_id = extract_id_from_header_line(line, 'FORMAT')
        if format_id and format_id not in existing_format_ids:
            merged['format_lines'].append(line)

    for line in octopus_defaults.get('alt_lines', []):
        alt_id = extract_id_from_header_line(line, 'ALT')
        if alt_id and alt_id not in existing_alt_ids:
            merged['alt_lines'].append(line)

    return merged


def generate_sv_header(contig_lines, input_vcf_file=None):
    """
    Generate SVCF file header lines according to SVCF specification.
    If input_vcf_file is provided, extract and preserve original header definitions.

    🔴 CHANGED: When the input VCF has more than one sample column, the
    generated header includes a ##OctopuSV_mode=multi marker and the #CHROM
    line lists all original sample names. Single-sample behavior is
    unchanged (the #CHROM line ends with the original sample name, falling
    back to "Sample" if no input file is provided).
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

    # If input VCF file is provided, extract original definitions
    if input_vcf_file:
        original_headers = extract_original_header_definitions(input_vcf_file)
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

    # 🔴 CHANGED: insert a multi-sample mode marker so downstream tools
    # (e.g. svcf2vcf) know to preserve all sample columns.
    final_header = list(basic_header)
    if is_multi_sample:
        final_header.append("##OctopuSV_mode=multi")

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


def write_sv_vcf(contig_lines, events, output_file, input_vcf_file=None):
    """
    Write SV events to VCF file with proper header definitions.
    If input_vcf_file is provided, preserve original header definitions.
    """
    sv_header = generate_sv_header(contig_lines, input_vcf_file)
    with open(output_file, "w") as f:
        for header in sv_header:
            f.write(header + "\n")
        for event in events:
            f.write(str(event) + "\n")