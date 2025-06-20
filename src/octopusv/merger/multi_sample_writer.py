import datetime
import os
from typing import List, Dict, Any
from .name_mapper import NameMapper


class MultiSampleWriter:
    """Handle writing VCF output for sample mode with consistent ordering."""

    def __init__(self, name_mapper: NameMapper):
        """Initialize writer with name mapper."""
        self.name_mapper = name_mapper

    def write_results(self, output_file, events, contigs, sv_merger):
        """Write merged results to output file with consistent SOURCES and SAMPLE ordering."""
        with open(output_file, "w") as f:
            self._write_header(f, contigs)

            for event in events:
                self._write_event(f, event, sv_merger)

    def _write_header(self, file_handle, contigs):
        """Write VCF header including sample columns in input file order."""
        # Write standard VCF header
        file_handle.write("##fileformat=SVCFv1.0\n")
        file_date = datetime.datetime.now().strftime("%Y-%m-%d|%I:%M:%S%p|")
        file_handle.write(f"##fileDate={file_date}\n")
        file_handle.write("##source=octopusV\n")

        # Write contig information
        for contig_id, contig_length in contigs.items():
            file_handle.write(f"##contig=<ID={contig_id},length={contig_length}>\n")

        # Write column headers for sample mode - ordered by input file sequence
        sample_names = self.name_mapper.get_all_display_names()
        header_line = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(sample_names) + "\n"
        file_handle.write(header_line)

    def _write_event(self, file_handle, event, sv_merger):
        """Write a single event record with consistent SOURCES and SAMPLE ordering."""
        # Step 1: Get ordered sources according to input file order
        ordered_sources = sv_merger._get_ordered_sources_for_event(event)

        # Step 2: Generate SOURCES field with consistent ordering
        display_sources = ",".join([self.name_mapper.get_display_name(f) for f in ordered_sources])

        # Step 3: Prepare INFO field with ordered SOURCES
        info_items = []
        for k, v in event.info.items():
            if k == "SOURCES":
                # Replace with our ordered sources using display names
                info_items.append(f"SOURCES={display_sources}")
            else:
                info_items.append(f"{k}={v}")

        info_field = ";".join(info_items)
        if "SOURCES" not in info_field:
            info_field += f";SOURCES={display_sources}"

        # Step 4: Get FORMAT field
        format_field = event.format
        format_keys = format_field.split(":")

        # Step 5: Format sample data for sample mode with consistent ordering
        sample_part = self._format_sample_mode_output(event, format_keys, sv_merger, ordered_sources)

        # Step 6: Write the complete record
        record_part1 = f"{event.chrom}\t{event.pos}\t{event.sv_id}\t{event.ref}\t{event.alt}\t"
        record_part2 = f"{event.quality}\t{event.filter}\t{info_field}\t{format_field}\t"
        file_handle.write(record_part1 + record_part2 + sample_part + "\n")

    def _format_sample_mode_output(self, event, format_keys, sv_merger, ordered_sources):
        """Format output for sample mode with consistent ordering based on input file sequence.

        Args:
            event: The SV event object
            format_keys: List of FORMAT field keys
            sv_merger: The SVMerger instance
            ordered_sources: List of source files in correct input order
        """
        sample_names = self.name_mapper.get_all_display_names()

        # Get caller mode output data ordered by input file sequence
        caller_mode_output = self._get_caller_mode_output_ordered(event, format_keys, sv_merger, ordered_sources)

        if not caller_mode_output or caller_mode_output == "./.":
            # No data available, fill with missing data for all samples
            missing_data = "0/0" + ":" + ":".join(["."] * (len(format_keys) - 1))
            return "\t".join([missing_data] * len(sample_names))

        # Split caller mode output into individual sample data pieces
        caller_data_parts = caller_mode_output.split('\t')

        # Map ordered sources to their corresponding display names and data
        source_data_map = {}
        for i, source_file in enumerate(ordered_sources):
            if i < len(caller_data_parts):
                sample_display_name = self.name_mapper.get_display_name(source_file)
                source_data_map[sample_display_name] = caller_data_parts[i]

        # Generate sample columns in the same order as header (input file order)
        sample_columns = []
        missing_data = "0/0" + ":" + ":".join(["."] * (len(format_keys) - 1))

        for sample_name in sample_names:
            if sample_name in source_data_map:
                sample_columns.append(source_data_map[sample_name])
            else:
                sample_columns.append(missing_data)

        return "\t".join(sample_columns)

    def _get_caller_mode_output_ordered(self, event, format_keys, sv_merger, ordered_sources):
        """Get caller mode output with data ordered according to input file sequence.

        Args:
            event: The SV event object
            format_keys: List of FORMAT field keys
            sv_merger: The SVMerger instance
            ordered_sources: List of source files in correct input order

        Returns:
            str: Tab-separated sample data ordered by input file sequence
        """
        merged_samples = getattr(event, "merged_samples", [])
        if merged_samples:
            # Reorder merged_samples to match ordered_sources
            reordered_samples = sv_merger._reorder_merged_samples(event, ordered_sources)

            # Format samples in the correct order
            sample_strings = []
            for _, _, sample_data in reordered_samples:
                if isinstance(sample_data, dict):
                    values = []
                    for key in format_keys:
                        value = sample_data.get(key, ".")
                        values.append(str(value))
                    sample_str = ":".join(values)
                    if sample_str.endswith(":.:."):
                        sample_str = sample_str[:-4]
                    sample_strings.append(sample_str)
                else:
                    sample_str = str(sample_data)
                    if sample_str.endswith(":.:."):
                        sample_str = sample_str[:-4]
                    sample_strings.append(sample_str)

            return "\t".join(sample_strings)
        elif hasattr(event, "sample"):
            # Single sample case
            formatted_values = sv_merger.format_sample_values(format_keys, event.sample)
            if formatted_values.endswith(":.:."):
                formatted_values = formatted_values[:-4]
            return formatted_values
        else:
            return "./."