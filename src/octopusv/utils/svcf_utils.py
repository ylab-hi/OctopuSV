from datetime import datetime


def generate_sv_header(contig_lines):
    """Generates SVCF file header lines according to SVCF specification, including original contig lines."""
    current_time_str = datetime.now().strftime("%Y-%m-%d|%I:%M:%S%p|%Z")
    return [
        "##fileformat=VCFv4.2",
        f"##fileDate={current_time_str}",  # Use current time
        "##source=OctopuSV",
        *contig_lines,  # Include original ##contig lines
        '##ALT=<ID=DEL,Description="Deletion">',
        '##ALT=<ID=INV,Description="Inversion">',
        '##ALT=<ID=INS,Description="Insertion">',
        '##ALT=<ID=DUP,Description="Duplication">',
        '##ALT=<ID=TRA,Description="Translocation">',
        '##ALT=<ID=BND,Description="Breakend">',
        '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">',
        '##INFO=<ID=CHR2,Number=1,Type=String,Description="Chromosome for end">',
        '##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant described in this record">',
        '##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Difference in length between REF and ALT alleles">',
        '##INFO=<ID=SUPPORT,Number=1,Type=Integer,Description="Number of pieces of evidence supporting the variant">',
        '##INFO=<ID=SVMETHOD,Number=1,Type=String,Description="The software used to identify the SV">',
        '##INFO=<ID=RTID,Number=1,Type=String,Description="Associated ID for reciprocal translocations if available">',
        '##INFO=<ID=AF,Number=1,Type=Float,Description="Allele Frequency">',
        '##INFO=<ID=STRAND,Number=1,Type=String,Description="Strand orientation of the SV">',
        '##FILTER=<ID=PASS,Description="All filters passed, variant is most likely true">',
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
        '##FORMAT=<ID=CO,Number=1,Type=String,Description="Coordinate information of the SV">',
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSample",
    ]


def write_sv_vcf(contig_lines, events, output_file):
    """Writes the SV events to a VCF file formatted according to the SVCF specification."""
    sv_header = generate_sv_header(contig_lines)  # Get the new SVCF headers, including contig lines

    with open(output_file, "w") as f:
        for header in sv_header:
            f.write(header + "\n")
        for event in events:
            f.write(str(event) + "\n")
