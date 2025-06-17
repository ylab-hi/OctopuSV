import datetime
import os
from typing import List, Dict, Any
from .name_mapper import NameMapper


class MultiSampleWriter:
    """Handle writing VCF output for sample mode."""

    def __init__(self, name_mapper: NameMapper):
        """Initialize writer with name mapper."""
        self.name_mapper = name_mapper

    def write_results(self, output_file, events, contigs, sv_merger):
        """Write merged results to output file."""
        with open(output_file, "w") as f:
            self._write_header(f, contigs)

            for event in events:
                self._write_event(f, event, sv_merger)

    def _write_header(self, file_handle, contigs):
        """Write VCF header including sample columns."""
        # Write standard VCF header
        file_handle.write("##fileformat=SVCFv1.0\n")
        file_date = datetime.datetime.now().strftime("%Y-%m-%d|%I:%M:%S%p|")
        file_handle.write(f"##fileDate={file_date}\n")
        file_handle.write("##source=octopusV\n")

        # Write contig information
        for contig_id, contig_length in contigs.items():
            file_handle.write(f"##contig=<ID={contig_id},length={contig_length}>\n")

        # Write column headers for sample mode
        sample_names = self.name_mapper.get_all_display_names()
        header_line = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(sample_names) + "\n"
        file_handle.write(header_line)

    def _write_event(self, file_handle, event, sv_merger):
        """Write a single event record."""
        # Prepare INFO field with converted source names
        info_items = []
        for k, v in event.info.items():
            if k == "SOURCES":
                display_sources = self.name_mapper.convert_source_string(v)
                info_items.append(f"SOURCES={display_sources}")
            else:
                info_items.append(f"{k}={v}")

        info_field = ";".join(info_items)
        if "SOURCES" not in info_field:
            display_sources = self.name_mapper.convert_source_string(event.source_file)
            info_field += f";SOURCES={display_sources}"

        # Get FORMAT field
        format_field = event.format
        format_keys = format_field.split(":")

        # Format sample data for sample mode
        sample_part = self._format_sample_mode_output(event, format_keys, sv_merger)

        # Write the record
        record_part1 = f"{event.chrom}\t{event.pos}\t{event.sv_id}\t{event.ref}\t{event.alt}\t"
        record_part2 = f"{event.quality}\t{event.filter}\t{info_field}\t{format_field}\t"
        file_handle.write(record_part1 + record_part2 + sample_part + "\n")

    def _format_sample_mode_output(self, event, format_keys, sv_merger):
        """Format output for sample mode - directly reuse caller mode logic."""
        sample_names = self.name_mapper.get_all_display_names()

        caller_mode_output = self._get_caller_mode_output(event, format_keys, sv_merger)

        if not caller_mode_output or caller_mode_output == "./.":
            missing_data = "0/0" + ":" + ":".join(["."] * (len(format_keys) - 1))
            return "\t".join([missing_data] * len(sample_names))

        caller_data_parts = caller_mode_output.split('\t')

        source_files = event.source_file.split(",")
        source_files = [f.strip() for f in source_files]

        sample_data_map = {}

        for i, source_file in enumerate(source_files):
            if i < len(caller_data_parts):
                sample_display_name = self.name_mapper.get_display_name(source_file)
                sample_data_map[sample_display_name] = caller_data_parts[i]

        sample_columns = []
        missing_data = "0/0" + ":" + ":".join(["."] * (len(format_keys) - 1))

        for sample_name in sample_names:
            if sample_name in sample_data_map:
                sample_columns.append(sample_data_map[sample_name])
            else:
                sample_columns.append(missing_data)

        return "\t".join(sample_columns)

    def _get_caller_mode_output(self, event, format_keys, sv_merger):
        """Get the exact same output that caller mode would produce."""
        merged_samples = getattr(event, "merged_samples", [])
        if merged_samples:
            # Find the representative sample and other samples
            rep_sample_data = None
            other_samples = []
            for sample_name, sample_format, sample_data in merged_samples:
                if isinstance(sample_data, dict) and sample_data.get("ID") == event.sv_id:
                    rep_sample_data = (sample_name, sample_format, sample_data)
                else:
                    other_samples.append((sample_name, sample_format, sample_data))

            # Format samples in the correct order (representative first, then others)
            sample_strings = []

            # Add representative sample first
            if rep_sample_data:
                _, _, sample_data = rep_sample_data
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
                    if str(sample_data).endswith(":.:."):
                        sample_data = str(sample_data)[:-4]
                    sample_strings.append(str(sample_data))

            # Add other samples
            for _, _, sample_data in other_samples:
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
                    if str(sample_data).endswith(":.:."):
                        sample_data = str(sample_data)[:-4]
                    sample_strings.append(str(sample_data))

            return "\t".join(sample_strings)
        elif hasattr(event, "sample"):
            formatted_values = sv_merger.format_sample_values(format_keys, event.sample)
            if formatted_values.endswith(":.:."):
                formatted_values = formatted_values[:-4]
            return formatted_values
        else:
            return "./."