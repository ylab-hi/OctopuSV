import logging
import re
import sys

from octopusv.sv import SVEvent


def is_same_chr_bnd(event):
    """Check if the POS and ALT of an event are on the same chromosome."""
    if event.is_BND():
        split_result = re.split(r"[\[\]:]", event.alt)
        if len(split_result) != 4:
            return False
        else:
            chrom_alt, _ = split_result[1:3]
            return event.chrom == chrom_alt

    return False  # For non-BND, we won't categorize them as same_chr_bnd or diff_chr_bnd events


def check_vcf_format(vcf_file_path):
    """Check the format of a VCF file.

    Allow both standard VCF (10 columns) and simplified VCF (8 columns) formats.
    Raise an error and exit if the format is incorrect.
    """
    with open(vcf_file_path) as f:
        lines = f.readlines()

    # Check if there is at least one header line
    if not any(line.startswith("#") for line in lines):
        logging.error(
            "Invalid VCF format. The file should contain at least one header line starting with '#'.",
        )
        sys.exit(1)

    for line in lines:
        if line.startswith("#"):
            continue  # Skip header lines

        # Check for space in lines
        if " " in line:
            logging.error(
                "Invalid VCF format. Non-header lines should not contain spaces.",
            )
            sys.exit(1)

        fields = line.strip().split("\t")

        # Check the minimum number of columns (8 for simplified VCF)
        if len(fields) < 8:
            logging.error(
                f"Invalid VCF format. Expected at least 8 fields, but got {len(fields)}",
            )
            sys.exit(1)

        # Check that the position is a number
        try:
            int(fields[1])
        except ValueError:
            logging.error(
                f"Invalid VCF format. Position (field 2) should be a number, but got {fields[1]}",
            )
            sys.exit(1)

        # Check that the quality score is a number or '.'
        if fields[5] != ".":
            try:
                float(fields[5])
            except ValueError:
                logging.error(
                    f"Invalid VCF format. Quality score (field 6) should be a number or '.', but got {fields[5]}",
                )
                sys.exit(1)


def parse_vcf(vcf_file_path):
    """Parse VCF file into lists of SVEvent objects based on their type.

    Handles standard VCF (10 columns), simplified VCF (8 columns) and
    multi-sample VCF (>10 columns) formats. For multi-sample VCFs each
    SVEvent retains all sample columns via its `samples` attribute, and
    they are propagated to the output SVCF.
    """
    check_vcf_format(vcf_file_path)  # Check the format first
    same_chr_bnd_events = []
    diff_chr_bnd_events = []
    non_bnd_events = []
    contig_lines = []  # Store ##contig lines here
    is_svaba_output = False  # Flag to detect if it's SVABA output
    source_info = "."  # Default value of source

    with open(vcf_file_path) as f:
        for line in f:
            if line.startswith("##source="):
                source_info = line.split("=")[1].split(" ")[0].strip()
                if "svaba" in line.lower():
                    is_svaba_output = True
            elif line.startswith("##contig"):
                contig_lines.append(line.strip())
            elif not line.startswith("#"):  # Skip all header lines except ##contig
                fields = line.strip().split("\t")

                # 🔴 CHANGED: Normalize fields into (core_fields, sample_fields).
                # core_fields = the first 9 VCF columns (CHROM..FORMAT).
                # sample_fields = list of one or more raw sample column strings.
                # This split is what lets us support multi-sample VCFs without
                # exploding the SVEvent positional-argument signature.
                if len(fields) == 8:
                    # Simplified VCF: synthesize a minimal FORMAT + single SAMPLE.
                    core_fields = fields[:8] + ["GT"]
                    sample_fields = ["0/1"]
                elif is_svaba_output and len(fields) == 13:
                    # SVABA special handling preserved: use original column 8
                    # (FORMAT) and column 12 as the single relevant sample.
                    core_fields = fields[:8] + [fields[8]]
                    sample_fields = [fields[12]]
                elif not is_svaba_output and len(fields) == 13:
                    raise ValueError("VCF format error: Detected 13 columns in the header.")
                elif len(fields) >= 10:
                    # Standard single-sample (len == 10) OR multi-sample (>10):
                    # split into 9 core columns + N sample columns.
                    core_fields = fields[:9]
                    sample_fields = fields[9:]
                elif len(fields) == 9:
                    # Has FORMAT but no SAMPLE column: pad with default '0/1'.
                    core_fields = fields[:9]
                    sample_fields = ["0/1"]
                else:
                    # Shouldn't reach here because check_vcf_format rejects <8.
                    continue

                # 🔴 CHANGED: pass `samples` so multi-sample info is preserved
                # all the way through the converter pipeline (the converters
                # use copy.deepcopy, so they carry `samples` automatically).
                event = SVEvent(
                    *core_fields,
                    sample=sample_fields[0],
                    samples=sample_fields,
                )
                event.source = source_info  # Add source dynamically

                # Skip non-variant records (e.g., Dragen REF regions without SVTYPE)
                if 'SVTYPE' not in event.info:
                    continue

                # Convert generic CNV to specific DEL/DUP based on ALT and ID fields
                if event.info.get('SVTYPE') == 'CNV':
                    # Method 1: Infer from ALT field
                    if '<DEL>' in event.alt.upper():
                        event.info['SVTYPE'] = 'DEL'
                    elif '<DUP>' in event.alt.upper():
                        event.info['SVTYPE'] = 'DUP'
                    # Method 2: Infer from ID field (Dragen uses LOSS/GAIN in ID)
                    elif 'LOSS' in event.id.upper():
                        event.info['SVTYPE'] = 'DEL'
                    elif 'GAIN' in event.id.upper():
                        event.info['SVTYPE'] = 'DUP'
                    # Method 3: Infer from SVLEN (negative=deletion, positive=duplication)
                    elif 'SVLEN' in event.info:
                        try:
                            svlen = int(event.info['SVLEN'])
                            if svlen < 0:
                                event.info['SVTYPE'] = 'DEL'
                            elif svlen > 0:
                                event.info['SVTYPE'] = 'DUP'
                        except ValueError:
                            pass
                    # If still CNV after all methods, log warning and skip
                    if event.info.get('SVTYPE') == 'CNV':
                        logging.warning(f"Could not determine specific type for CNV event {event.id}, skipping")
                        continue

                if event.is_BND():
                    if is_same_chr_bnd(event):  # Check if the event is same chromosome
                        same_chr_bnd_events.append(event)
                    else:
                        diff_chr_bnd_events.append(event)  # Different chromosome
                else:
                    non_bnd_events.append(event)  # Non-BND events

    return contig_lines, same_chr_bnd_events, diff_chr_bnd_events, non_bnd_events