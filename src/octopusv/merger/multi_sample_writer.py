import datetime
import os
from typing import List, Dict, Any
from .name_mapper import NameMapper


class MultiSampleWriter:
    """Handle writing VCF output for sample mode with consistent ordering."""

    def __init__(self, name_mapper: NameMapper):
        """Initialize writer with name mapper."""
        self.name_mapper = name_mapper

    def write_results(self, output_file, events, contigs):
        """Write merged results to output file with consistent SOURCES and SAMPLE ordering.

        Args:
            output_file: Path to output file
            events: List of preprocessed events (with ordered_samples attribute)
            contigs: Dictionary of contig information
        """
        with open(output_file, "w") as f:
            self._write_header(f, contigs)

            for event in events:
                self._write_event(f, event)

    def _write_header(self, file_handle, contigs):
        """Write VCF header including sample columns in input file order."""
        # Write standard VCF header
        file_handle.write("##fileformat=VCFv4.2\n")

        # Mark this as multi-sample mode for downstream tools
        file_handle.write("##OctopuSV_mode=multi\n")

        file_date = datetime.datetime.now().strftime("%Y-%m-%d|%I:%M:%S%p|")
        file_handle.write(f"##fileDate={file_date}\n")
        file_handle.write("##source=OctopuSV\n")
        file_handle.write(
            "##OctopuSV_WARNING=This is SVCF format. Use 'octopusv svcf2vcf' to change back to standard VCF format before bcftools/vcftools\n")

        # Write contig information
        for contig_id, contig_length in contigs.items():
            file_handle.write(f"##contig=<ID={contig_id},length={contig_length}>\n")

        # Write standard INFO/FORMAT definitions (minimal set)
        self._write_standard_definitions(file_handle)

        # Write column headers for sample mode - ordered by input file sequence
        sample_names = self.name_mapper.get_all_display_names()
        header_line = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(sample_names) + "\n"
        file_handle.write(header_line)

    def _write_standard_definitions(self, file_handle):
        """Write standard INFO and FORMAT field definitions."""
        # INFO definitions
        file_handle.write('##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">\n')
        file_handle.write('##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant">\n')
        file_handle.write(
            '##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Difference in length between REF and ALT alleles">\n')
        file_handle.write('##INFO=<ID=CHR2,Number=1,Type=String,Description="Chromosome for end coordinate">\n')
        file_handle.write(
            '##INFO=<ID=SUPPORT,Number=1,Type=Integer,Description="Number of reads supporting this variant">\n')
        file_handle.write('##INFO=<ID=SVMETHOD,Number=1,Type=String,Description="Method used to detect SV">\n')
        file_handle.write('##INFO=<ID=STRAND,Number=1,Type=String,Description="Strand orientation of the SV">\n')
        file_handle.write('##INFO=<ID=RNAMES,Number=.,Type=String,Description="Supporting read names">\n')
        file_handle.write(
            '##INFO=<ID=SOURCES,Number=.,Type=String,Description="List of input samples/files that support this variant">\n')

        # FORMAT definitions
        file_handle.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
        file_handle.write(
            '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths for the ref and alt alleles">\n')
        file_handle.write('##FORMAT=<ID=LN,Number=1,Type=Integer,Description="Length of SV">\n')
        file_handle.write('##FORMAT=<ID=ST,Number=1,Type=String,Description="Strand orientation of SV">\n')
        file_handle.write('##FORMAT=<ID=QV,Number=1,Type=Integer,Description="Quality value">\n')
        file_handle.write('##FORMAT=<ID=TY,Number=1,Type=String,Description="Type of SV">\n')
        file_handle.write('##FORMAT=<ID=ID,Number=1,Type=String,Description="Unique identifier for the SV">\n')
        file_handle.write('##FORMAT=<ID=SC,Number=1,Type=String,Description="Source caller/method">\n')
        file_handle.write('##FORMAT=<ID=REF,Number=1,Type=String,Description="Reference allele sequence">\n')
        file_handle.write('##FORMAT=<ID=ALT,Number=1,Type=String,Description="Alternate allele sequence">\n')
        file_handle.write('##FORMAT=<ID=CO,Number=1,Type=String,Description="Coordinate information of the SV">\n')

    def _write_event(self, file_handle, event):
        """Write a single event record with consistent SOURCES and SAMPLE ordering.

        Args:
            file_handle: File handle to write to
            event: Event object with ordered_samples attribute
        """
        # Get ordered sample names for SOURCES field
        all_sample_names = self.name_mapper.get_all_display_names()
        ordered_samples = getattr(event, "ordered_samples", [])

        # Determine which samples have this SV for SOURCES field
        sources = []
        for i, sample_data in enumerate(ordered_samples):
            if sample_data is not None:
                sources.append(all_sample_names[i])

        sources_str = ",".join(sources) if sources else "."

        # Prepare INFO field
        info_items = []
        for k, v in event.info.items():
            if k != "SOURCES":  # We'll add SOURCES explicitly
                info_items.append(f"{k}={v}")

        info_items.append(f"SOURCES={sources_str}")
        info_field = ";".join(info_items)

        # Get FORMAT field
        format_field = event.format
        format_keys = format_field.split(":")

        # Format sample data using preprocessed ordered_samples
        sample_part = self._format_sample_columns(ordered_samples, format_keys)

        # Write the complete record
        record = (f"{event.chrom}\t{event.pos}\t{event.sv_id}\t{event.ref}\t{event.alt}\t"
                  f"{event.quality}\t{event.filter}\t{info_field}\t{format_field}\t{sample_part}\n")
        file_handle.write(record)

    def _format_sample_columns(self, ordered_samples, format_keys):
        """Format sample columns using preprocessed ordered_samples.

        Args:
            ordered_samples: List of sample data (or None) in input file order
            format_keys: List of FORMAT field keys

        Returns:
            str: Tab-separated sample column data
        """
        sample_columns = []

        # Missing data default value (reference genotype)
        missing_data = "0/0" + ":" + ":".join(["."] * (len(format_keys) - 1))

        for sample_data in ordered_samples:
            if sample_data is not None:
                if isinstance(sample_data, dict):
                    values = [str(sample_data.get(key, ".")) for key in format_keys]
                    sample_str = ":".join(values)
                    # Clean up trailing :.:.
                    if sample_str.endswith(":.:."):
                        sample_str = sample_str[:-4]
                    sample_columns.append(sample_str)
                else:
                    sample_str = str(sample_data)
                    if sample_str.endswith(":.:."):
                        sample_str = sample_str[:-4]
                    sample_columns.append(sample_str)
            else:
                # This sample doesn't have this SV
                sample_columns.append(missing_data)

        return "\t".join(sample_columns)