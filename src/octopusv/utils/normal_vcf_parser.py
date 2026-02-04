import logging
import re
import sys

from octopusv.sv import SVEvent


def is_same_chr_bnd(event):
    """Check if the POS and ALT of an event are on the same chromosome."""
    if event.is_BND():
        split_result = re.split(r"[\[\]:]", event.alt)
        if len(split_result) != 4:
            logging.info(
                f"Unexpected ALT format, it should be something like N]chr10:69650962]: {split_result}",
            )
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

    Handles both standard VCF (10 columns) and simplified VCF (8 columns) formats.
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

                # Add default FORMAT and SAMPLE fields if they don't exist
                if len(fields) == 8:
                    fields.extend(["GT", "0/1"])  # Add minimal FORMAT and SAMPLE fields

                # Adjust for SVABA if it's detected by the source line and has 13 fields
                if is_svaba_output and len(fields) == 13:
                    adjusted_fields = fields[:8] + [fields[8], fields[12]]
                elif not is_svaba_output and len(fields) == 13:
                    raise ValueError("VCF format error: Detected 13 columns in the header.")
                else:
                    adjusted_fields = fields

                event = SVEvent(*adjusted_fields)  # Unpack fields and send to SVEvent class
                event.source = source_info  # Add source dynamically

                if event.is_BND():
                    if is_same_chr_bnd(event):  # Check if the event is same chromosome
                        same_chr_bnd_events.append(event)
                    else:
                        diff_chr_bnd_events.append(event)  # Different chromosome
                else:
                    non_bnd_events.append(event)  # Non-BND events

    return contig_lines, same_chr_bnd_events, diff_chr_bnd_events, non_bnd_events
