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
        file_handle.write("##fileformat=VCFv4.2\n")
        file_date = datetime.datetime.now().strftime("%Y-%m-%d|%I:%M:%S%p|")
        file_handle.write(f"##fileDate={file_date}\n")
        file_handle.write("##source=OctopuSV\n")

        # Write contig information
        for contig_id, contig_length in contigs.items():
            file_handle.write(f"##contig=<ID={contig_id},length={contig_length}>\n")

        # Write column headers for sample mode - ordered by input file sequence
        sample_names = self.name_mapper.get_all_display_names()
        header_line = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(sample_names) + "\n"
        file_handle.write(header_line)

    def _write_event(self, file_handle, event, sv_merger):
        """Write a single event record with consistent SOURCES and SAMPLE ordering."""
        # Step 1: Get ordered sources according to input file order (internal method)
        ordered_sources = self._get_ordered_sources_for_event(event, sv_merger.all_input_files)

        # Step 2: Generate SOURCES field with consistent ordering
        display_sources = ",".join([self.name_mapper.get_display_name(f) for f in ordered_sources])

        # Step 3: Prepare INFO field with ordered SOURCES
        info_items = []
        for k, v in event.info.items():
            if k == "SOURCES":
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

    def _get_ordered_sources_for_event(self, event, all_input_files):
        """Get source files for an event ordered according to input file order.
        Internal method - doesn't require changes to SVMerger.

        Args:
            event: The SV event object
            all_input_files: List of all input files from SVMerger

        Returns:
            list: List of source file paths ordered by input file order
        """
        # Get source files from the event
        source_files = event.source_file.split(",")
        source_files = [f.strip() for f in source_files]

        # Order them according to the input file order
        ordered_sources = []
        for input_file in all_input_files:
            input_basename = os.path.basename(input_file)
            # Check if this input file is among the event's sources
            for source_file in source_files:
                source_basename = os.path.basename(source_file)
                if input_basename == source_basename:
                    ordered_sources.append(source_file)
                    break

        return ordered_sources

    def _format_sample_mode_output(self, event, format_keys, sv_merger, ordered_sources):
        """Format output for sample mode with consistent ordering based on input file sequence.

        Args:
            event: The SV event object
            format_keys: List of FORMAT field keys
            sv_merger: The SVMerger instance
            ordered_sources: List of source files in correct input order
        """
        all_sample_names = self.name_mapper.get_all_display_names()
        merged_samples = getattr(event, "merged_samples", [])

        # Create a mapping from source basename to sample data
        source_to_sample = {}
        for sample_name, sample_format, sample_data in merged_samples:
            # Try to determine which source file this sample came from
            sample_assigned = False
            for input_file in sv_merger.all_input_files:
                input_basename = os.path.basename(input_file)

                # Check multiple ways to match sample to source
                if isinstance(sample_data, dict):
                    # Method 1: Check source_file field
                    source_file_field = sample_data.get('source_file', '')
                    if input_basename in source_file_field or input_file in source_file_field:
                        source_to_sample[input_basename] = sample_data
                        sample_assigned = True
                        break

                    # Method 2: Check if input file appears in sample data
                    sample_str = str(sample_data)
                    input_name = os.path.splitext(input_basename)[0]
                    if input_name.lower() in sample_str.lower():
                        source_to_sample[input_basename] = sample_data
                        sample_assigned = True
                        break

            # If no specific match found, assign to first available input file that contributed to this event
            if not sample_assigned:
                event_source_basenames = {os.path.basename(f.strip()) for f in event.source_file.split(",")}
                for input_file in sv_merger.all_input_files:
                    input_basename = os.path.basename(input_file)
                    if input_basename not in source_to_sample and input_basename in event_source_basenames:
                        source_to_sample[input_basename] = sample_data
                        break

        # Generate sample columns in the same order as header (input file order)
        sample_columns = []
        missing_data = "0/0" + ":" + ":".join(["."] * (len(format_keys) - 1))

        for i, sample_name in enumerate(all_sample_names):
            # Get corresponding input file for this sample position
            if i < len(sv_merger.all_input_files):
                input_file = sv_merger.all_input_files[i]
                input_basename = os.path.basename(input_file)

                if input_basename in source_to_sample:
                    sample_data = source_to_sample[input_basename]
                    if isinstance(sample_data, dict):
                        values = []
                        for key in format_keys:
                            value = sample_data.get(key, ".")
                            values.append(str(value))
                        sample_str = ":".join(values)
                        if sample_str.endswith(":.:."):
                            sample_str = sample_str[:-4]
                        sample_columns.append(sample_str)
                    else:
                        sample_str = str(sample_data)
                        if sample_str.endswith(":.:."):
                            sample_str = sample_str[:-4]
                        sample_columns.append(sample_str)
                else:
                    # This sample doesn't have this SV, use 0/0 (reference genotype)
                    sample_columns.append(missing_data)
            else:
                sample_columns.append(missing_data)

        return "\t".join(sample_columns)