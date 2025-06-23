import datetime
import os
import re

from .sv_merge_logic import should_merge
from .sv_selector import select_representative_sv
from .tra_merger import TRAMerger
from .bnd_merger import BNDMerger


class SVMerger:
    def __init__(
            self,
            classified_events,
            all_input_files,
            tra_delta=50,
            tra_min_overlap_ratio=0.5,
            tra_strand_consistency=True,
            max_distance=50,
            max_length_ratio=1.3,
            min_jaccard=0.7,
            bnd_delta=50,
    ):
        """Initialize SVMerger with the given parameters and events.

        Args:
            classified_events: Dictionary of classified SV events
            all_input_files: List of all input file paths
            tra_delta: Position uncertainty threshold for TRA events
            tra_min_overlap_ratio: Minimum overlap ratio for TRA events
            tra_strand_consistency: Whether to require strand consistency for TRA events
            max_distance: Maximum allowed distance between start or end positions
            max_length_ratio: Maximum allowed ratio between event lengths
            min_jaccard: Minimum required Jaccard index for overlap
            bnd_delta: Position uncertainty threshold for BND events (default: 50)
        """
        self.classified_events = classified_events
        self.all_input_files = [str(file) for file in all_input_files]
        self.merged_events: dict[str, dict[str, list]] = {}
        self.event_groups: dict[str, dict[str, list[list]]] = {}
        self.tra_merger = TRAMerger(tra_delta, tra_min_overlap_ratio, tra_strand_consistency)
        self.bnd_merger = BNDMerger(bnd_delta)
        self.max_distance = max_distance
        self.max_length_ratio = max_length_ratio
        self.min_jaccard = min_jaccard

    def merge(self):
        """Merge all SV events based on their types and chromosomes."""
        for sv_type, chromosomes in self.classified_events.items():
            if sv_type == "TRA":
                for (_chr1, _chr2), events in chromosomes.items():
                    for event in events:
                        self.tra_merger.add_event(event)
            elif sv_type == "BND":
                for chromosome, events in chromosomes.items():
                    for event in events:
                        self.bnd_merger.add_event(event)
            else:
                if sv_type not in self.merged_events:
                    self.merged_events[sv_type] = {}
                    self.event_groups[sv_type] = {}
                for chromosome, events in chromosomes.items():
                    if chromosome not in self.merged_events[sv_type]:
                        self.merged_events[sv_type][chromosome] = []
                        self.event_groups[sv_type][chromosome] = []
                    for event in sorted(events, key=lambda e: (e.start_pos, e.end_pos, e.sv_id)):
                        self.add_and_merge_event(sv_type, chromosome, event)

    def add_and_merge_event(self, sv_type, chromosome, new_event):
        """Add a new event and merge it with existing events if possible."""
        events = self.merged_events[sv_type][chromosome]
        event_groups = self.event_groups[sv_type][chromosome]
        for idx, existing_event in enumerate(events):
            if should_merge(existing_event, new_event, self.max_distance, self.max_length_ratio, self.min_jaccard):
                event_groups[idx].append(new_event)
                return
        events.append(new_event)
        event_groups.append([new_event])

    def get_events(self, sv_type, chromosome, start, end):
        """Get events of given type within the specified region."""
        if sv_type == "TRA":
            return self.tra_merger.get_merged_events()
        elif sv_type == "BND":
            return self.bnd_merger.get_merged_events()
        if sv_type in self.event_groups and chromosome in self.event_groups[sv_type]:
            events = []
            for sv_group in self.event_groups[sv_type][chromosome]:
                representative_sv = select_representative_sv(sv_group)
                if representative_sv.start_pos <= end and representative_sv.end_pos >= start:
                    events.append(representative_sv)
            return events
        return []

    def get_all_merged_events(self):
        """Get all merged events across all types and chromosomes."""
        merged_events = []
        merged_events.extend(self.tra_merger.get_merged_events())
        merged_events.extend(self.bnd_merger.get_merged_events())
        for sv_type in self.event_groups:
            for chromosome in self.event_groups[sv_type]:
                for sv_group in self.event_groups[sv_type][chromosome]:
                    representative_sv = select_representative_sv(sv_group)
                    merged_events.append(representative_sv)
        return merged_events

    def get_events_by_source(self, sources, operation="union"):
        """Get events based on their source files and specified operation."""
        tra_events = self.tra_merger.get_merged_events()
        bnd_events = self.bnd_merger.get_merged_events()
        other_events = self.get_all_merged_events()
        other_events = [e for e in other_events if e.sv_type not in ["TRA", "BND"]]
        sources_set = {os.path.basename(source) for source in sources}

        if operation == "union":
            tra_filtered = [
                event for event in tra_events
                if sources_set.intersection({os.path.basename(s) for s in event.source_file.split(",")})
            ]
            bnd_filtered = [
                event for event in bnd_events
                if sources_set.intersection({os.path.basename(s) for s in event.source_file.split(",")})
            ]
            other_filtered = [
                event for event in other_events
                if sources_set.intersection({os.path.basename(s) for s in event.source_file.split(",")})
            ]
        elif operation == "intersection":
            tra_filtered = [
                event for event in tra_events
                if sources_set.issubset({os.path.basename(s) for s in event.source_file.split(",")})
            ]
            bnd_filtered = [
                event for event in bnd_events
                if sources_set.issubset({os.path.basename(s) for s in event.source_file.split(",")})
            ]
            other_filtered = [
                event for event in other_events
                if sources_set.issubset({os.path.basename(s) for s in event.source_file.split(",")})
            ]
        elif operation == "specific":
            source_file = next(iter(sources_set))
            other_files = {os.path.basename(f) for f in self.all_input_files} - {source_file}
            tra_filtered = [
                event for event in tra_events
                if source_file in [os.path.basename(s) for s in event.source_file.split(",")] and not any(
                    other in [os.path.basename(s) for s in event.source_file.split(",")] for other in other_files
                )
            ]
            bnd_filtered = [
                event for event in bnd_events
                if source_file in [os.path.basename(s) for s in event.source_file.split(",")] and not any(
                    other in [os.path.basename(s) for s in event.source_file.split(",")] for other in other_files
                )
            ]
            other_filtered = [
                event for event in other_events
                if source_file in [os.path.basename(s) for s in event.source_file.split(",")] and not any(
                    other in [os.path.basename(s) for s in event.source_file.split(",")] for other in other_files
                )
            ]
        else:
            msg = f"Unsupported operation: {operation}"
            raise ValueError(msg)

        return other_filtered + tra_filtered + bnd_filtered

    def get_events_by_exact_support(self, exact_support):
        """Get events supported by exactly N files."""
        tra_events = self.tra_merger.get_merged_events()
        bnd_events = self.bnd_merger.get_merged_events()
        other_events = self.get_all_merged_events()
        other_events = [e for e in other_events if e.sv_type not in ["TRA", "BND"]]

        tra_filtered = [event for event in tra_events if len(set(event.source_file.split(","))) == exact_support]
        bnd_filtered = [event for event in bnd_events if len(set(event.source_file.split(","))) == exact_support]
        other_filtered = [event for event in other_events if len(set(event.source_file.split(","))) == exact_support]
        return other_filtered + tra_filtered + bnd_filtered

    def get_events_by_support_range(self, min_support=None, max_support=None):
        """Get events supported by a range of files."""
        tra_events = self.tra_merger.get_merged_events()
        bnd_events = self.bnd_merger.get_merged_events()
        other_events = self.get_all_merged_events()
        other_events = [e for e in other_events if e.sv_type not in ["TRA", "BND"]]

        def within_range(event):
            support_count = len(set(event.source_file.split(",")))
            if min_support is not None and support_count < min_support:
                return False
            return not (max_support is not None and support_count > max_support)

        tra_filtered = [event for event in tra_events if within_range(event)]
        bnd_filtered = [event for event in bnd_events if within_range(event)]
        other_filtered = [event for event in other_events if within_range(event)]
        return other_filtered + tra_filtered + bnd_filtered

    def get_events_by_expression(self, expression):
        """Get events that satisfy a logical expression."""
        tra_events = self.tra_merger.get_merged_events()
        bnd_events = self.bnd_merger.get_merged_events()
        other_events = self.get_all_merged_events()
        other_events = [e for e in other_events if e.sv_type not in ["TRA", "BND"]]

        tra_filtered = [
            event for event in tra_events
            if self.evaluate_expression(expression, [os.path.basename(s) for s in event.source_file.split(",")])
        ]
        bnd_filtered = [
            event for event in bnd_events
            if self.evaluate_expression(expression, [os.path.basename(s) for s in event.source_file.split(",")])
        ]
        other_filtered = [
            event for event in other_events
            if self.evaluate_expression(expression, [os.path.basename(s) for s in event.source_file.split(",")])
        ]
        return other_filtered + tra_filtered + bnd_filtered

    def evaluate_expression(self, expression, event_sources):
        """Evaluate a logical expression against event sources."""

        def make_identifier(file_path):
            file_name = os.path.basename(file_path)
            return re.sub(r"\W|^(?=\d)", "_", file_name)

        context = {}
        all_sources = set(self.all_input_files)
        for source in all_sources:
            identifier = make_identifier(source)
            context[identifier] = False

        for source in event_sources:
            identifier = make_identifier(source)
            context[identifier] = True

        expr = expression
        for source in all_sources:
            identifier = make_identifier(source)
            expr = re.sub(r"\b" + re.escape(os.path.basename(source)) + r"\b", identifier, expr)

        expr = expr.replace("AND", "and").replace("OR", "or").replace("NOT", "not")

        try:
            return eval(expr, {"__builtins__": {}}, context)
        except Exception as e:
            msg = f"Invalid expression: {e}"
            raise ValueError(msg)

    def format_sample_values(self, format_keys, sample_dict):
        """Format sample values according to the FORMAT field."""
        values = []
        for key in format_keys:
            value = sample_dict.get(key, ".")
            if isinstance(value, list):
                value = ",".join(map(str, value))
            elif value is None:
                value = "."
            values.append(str(value))

        result = ":".join(values)
        if result.endswith(":.:."):
            result = result[:-4]
        return result

    def write_results(self, output_file, events, contigs, mode="caller", name_mapper=None):
        """Write merged results to output file with SOURCES and SAMPLE ordering by input file order."""
        if name_mapper and mode == "sample":
            # Use multi-sample writer for sample mode
            from .multi_sample_writer import MultiSampleWriter
            writer = MultiSampleWriter(name_mapper)
            writer.write_results(output_file, events, contigs, self)
        else:
            # Enhanced caller mode with simple input file order-based logic
            with open(output_file, "w") as f:
                # Write VCF header
                f.write("##fileformat=VCFv4.2\n")
                file_date = datetime.datetime.now().strftime("%Y-%m-%d|%I:%M:%S%p|")
                f.write(f"##fileDate={file_date}\n")
                f.write("##source=OctopuSV\n")

                # Write contig information
                for contig_id, contig_length in contigs.items():
                    f.write(f"##contig=<ID={contig_id},length={contig_length}>\n")

                # Write column headers
                sample_names = ["SAMPLE"]
                header_line = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(sample_names) + "\n"
                f.write(header_line)

                for event in events:
                    # Step 1: Find which input files contributed to this event
                    event_source_files = event.source_file.split(",")
                    event_source_basenames = {os.path.basename(f.strip()) for f in event_source_files}

                    # Step 2: Build SOURCES and samples in input file order
                    sources_in_order = []
                    samples_in_order = []
                    merged_samples = getattr(event, "merged_samples", [])

                    # Create a mapping from source basename to sample data
                    source_to_sample = {}
                    for sample_name, sample_format, sample_data in merged_samples:
                        # Try to determine which source file this sample came from
                        sample_assigned = False
                        for input_file in self.all_input_files:
                            input_basename = os.path.basename(input_file)

                            # Check multiple ways to match sample to source
                            if isinstance(sample_data, dict):
                                # Method 1: Check source_file field
                                source_file_field = sample_data.get('source_file', '')
                                if input_basename in source_file_field or input_file in source_file_field:
                                    source_to_sample[input_basename] = (sample_name, sample_format, sample_data)
                                    sample_assigned = True
                                    break

                                # Method 2: Check if input file appears in sample data
                                sample_str = str(sample_data)
                                input_name = os.path.splitext(input_basename)[0]
                                if input_name.lower() in sample_str.lower():
                                    source_to_sample[input_basename] = (sample_name, sample_format, sample_data)
                                    sample_assigned = True
                                    break

                        # If no specific match found, assign to first available input file
                        if not sample_assigned:
                            for input_file in self.all_input_files:
                                input_basename = os.path.basename(input_file)
                                if input_basename not in source_to_sample and input_basename in event_source_basenames:
                                    source_to_sample[input_basename] = (sample_name, sample_format, sample_data)
                                    break

                    # Step 3: Generate SOURCES and samples in input file order
                    for input_file in self.all_input_files:
                        input_basename = os.path.basename(input_file)
                        if input_basename in event_source_basenames:
                            # This input file contributed to this event
                            if name_mapper:
                                display_name = name_mapper.get_display_name(input_file)
                            else:
                                display_name = os.path.splitext(input_basename)[0]
                            sources_in_order.append(display_name)

                            # Get corresponding sample data
                            if input_basename in source_to_sample:
                                samples_in_order.append(source_to_sample[input_basename])

                    # Step 4: Generate SOURCES field
                    display_sources = ",".join(sources_in_order)

                    # Step 5: Prepare INFO field with ordered SOURCES
                    info_items = []
                    for k, v in event.info.items():
                        if k == "SOURCES":
                            info_items.append(f"SOURCES={display_sources}")
                        else:
                            info_items.append(f"{k}={v}")

                    info_field = ";".join(info_items)
                    if "SOURCES" not in info_field:
                        info_field += f";SOURCES={display_sources}"

                    # Step 6: Get FORMAT field
                    format_field = event.format
                    format_keys = format_field.split(":")

                    # Step 7: Generate sample data in input file order
                    if samples_in_order:
                        sample_strings = []
                        for _, _, sample_data in samples_in_order:
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

                        sample_part = "\t".join(sample_strings)
                    elif hasattr(event, "sample"):
                        # Single sample case
                        formatted_values = self.format_sample_values(format_keys, event.sample)
                        if formatted_values.endswith(":.:."):
                            formatted_values = formatted_values[:-4]
                        sample_part = formatted_values
                    else:
                        sample_part = "./."

                    # Step 8: Write the complete record
                    record_part1 = f"{event.chrom}\t{event.pos}\t{event.sv_id}\t{event.ref}\t{event.alt}\t"
                    record_part2 = f"{event.quality}\t{event.filter}\t{info_field}\t{format_field}\t"
                    f.write(record_part1 + record_part2 + sample_part + "\n")