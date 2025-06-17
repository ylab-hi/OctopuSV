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
            bnd_delta=50,  # Add BND-specific delta parameter
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
        self.bnd_merger = BNDMerger(bnd_delta)  # Add BND merger
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
                # BND events are classified by single chromosome, not chromosome pairs
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
                    for event in events:
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

        # Add TRA events
        merged_events.extend(self.tra_merger.get_merged_events())

        # Add BND events
        merged_events.extend(self.bnd_merger.get_merged_events())

        # Add other SV types (DEL, DUP, INV, INS)
        for sv_type in self.event_groups:
            for chromosome in self.event_groups[sv_type]:
                for sv_group in self.event_groups[sv_type][chromosome]:
                    representative_sv = select_representative_sv(sv_group)
                    merged_events.append(representative_sv)
        return merged_events

    def get_events_by_source(self, sources, operation="union"):
        """Get events based on their source files and specified operation.

        Args:
            sources: List of source file paths
            operation: One of "union", "intersection", or "specific"
        """
        tra_events = self.tra_merger.get_merged_events()
        bnd_events = self.bnd_merger.get_merged_events()  # Add BND events
        other_events = self.get_all_merged_events()

        # Remove TRA and BND events from other_events to avoid duplication
        other_events = [e for e in other_events if e.sv_type not in ["TRA", "BND"]]

        # Normalize sources to basenames
        sources_set = {os.path.basename(source) for source in sources}

        if operation == "union":
            tra_filtered = [
                event
                for event in tra_events
                if sources_set.intersection({os.path.basename(s) for s in event.source_file.split(",")})
            ]
            bnd_filtered = [
                event
                for event in bnd_events
                if sources_set.intersection({os.path.basename(s) for s in event.source_file.split(",")})
            ]
            other_filtered = [
                event
                for event in other_events
                if sources_set.intersection({os.path.basename(s) for s in event.source_file.split(",")})
            ]
        elif operation == "intersection":
            tra_filtered = [
                event
                for event in tra_events
                if sources_set.issubset({os.path.basename(s) for s in event.source_file.split(",")})
            ]
            bnd_filtered = [
                event
                for event in bnd_events
                if sources_set.issubset({os.path.basename(s) for s in event.source_file.split(",")})
            ]
            other_filtered = [
                event
                for event in other_events
                if sources_set.issubset({os.path.basename(s) for s in event.source_file.split(",")})
            ]
        elif operation == "specific":
            # Get the target file and all other files
            source_file = next(iter(sources_set))
            other_files = {os.path.basename(f) for f in self.all_input_files} - {source_file}

            tra_filtered = [
                event
                for event in tra_events
                if source_file in [os.path.basename(s) for s in event.source_file.split(",")] and not any(
                    other in [os.path.basename(s) for s in event.source_file.split(",")] for other in other_files
                )
            ]
            bnd_filtered = [
                event
                for event in bnd_events
                if source_file in [os.path.basename(s) for s in event.source_file.split(",")] and not any(
                    other in [os.path.basename(s) for s in event.source_file.split(",")] for other in other_files
                )
            ]
            other_filtered = [
                event
                for event in other_events
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
        bnd_events = self.bnd_merger.get_merged_events()  # Add BND events
        other_events = self.get_all_merged_events()

        # Remove TRA and BND events from other_events to avoid duplication
        other_events = [e for e in other_events if e.sv_type not in ["TRA", "BND"]]

        tra_filtered = [event for event in tra_events if len(set(event.source_file.split(","))) == exact_support]
        bnd_filtered = [event for event in bnd_events if len(set(event.source_file.split(","))) == exact_support]
        other_filtered = [event for event in other_events if len(set(event.source_file.split(","))) == exact_support]
        return other_filtered + tra_filtered + bnd_filtered

    def get_events_by_support_range(self, min_support=None, max_support=None):
        """Get events supported by a range of files."""
        tra_events = self.tra_merger.get_merged_events()
        bnd_events = self.bnd_merger.get_merged_events()  # Add BND events
        other_events = self.get_all_merged_events()

        # Remove TRA and BND events from other_events to avoid duplication
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
        bnd_events = self.bnd_merger.get_merged_events()  # Add BND events
        other_events = self.get_all_merged_events()

        # Remove TRA and BND events from other_events to avoid duplication
        other_events = [e for e in other_events if e.sv_type not in ["TRA", "BND"]]

        tra_filtered = [
            event
            for event in tra_events
            if self.evaluate_expression(expression, [os.path.basename(s) for s in event.source_file.split(",")])
        ]
        bnd_filtered = [
            event
            for event in bnd_events
            if self.evaluate_expression(expression, [os.path.basename(s) for s in event.source_file.split(",")])
        ]
        other_filtered = [
            event
            for event in other_events
            if self.evaluate_expression(expression, [os.path.basename(s) for s in event.source_file.split(",")])
        ]
        return other_filtered + tra_filtered + bnd_filtered

    def evaluate_expression(self, expression, event_sources):
        """Evaluate a logical expression against event sources."""

        # Function to convert file names to valid Python identifiers
        def make_identifier(file_path):
            # Extract basename to ensure consistency
            file_name = os.path.basename(file_path)
            # Replace non-alphanumeric characters with underscores
            return re.sub(r"\W|^(?=\d)", "_", file_name)

        # Build context for all input files
        context = {}
        all_sources = set(self.all_input_files)
        for source in all_sources:
            identifier = make_identifier(source)
            context[identifier] = False  # Initialize all to False

        # Set True for event_sources
        for source in event_sources:
            identifier = make_identifier(source)
            context[identifier] = True

        # Replace file names in expression with identifiers
        expr = expression
        for source in all_sources:
            identifier = make_identifier(source)
            # Use regex to replace exact file names to avoid partial matches
            expr = re.sub(r"\b" + re.escape(os.path.basename(source)) + r"\b", identifier, expr)

        # Convert logical operators to Python syntax
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

        # Remove the :.: at the end if it exists
        result = ":".join(values)
        if result.endswith(":.:."):
            result = result[:-4]
        return result

    def write_results(self, output_file, events, contigs, mode="caller", name_mapper=None):
        """Write merged results to output file."""
        if name_mapper and mode == "sample":
            # Use new multi-sample writer
            from .multi_sample_writer import MultiSampleWriter
            writer = MultiSampleWriter(name_mapper)
            writer.write_results(output_file, events, contigs, self)
        else:
            # Use existing logic for caller mode (keep all current code unchanged)
            with open(output_file, "w") as f:
                # Write VCF header
                f.write("##fileformat=VCFv4.2\n")
                file_date = datetime.datetime.now().strftime("%Y-%m-%d|%I:%M:%S%p|")
                f.write(f"##fileDate={file_date}\n")
                f.write("##source=OctopuSV\n")

                # Write contig information
                for contig_id, contig_length in contigs.items():
                    f.write(f"##contig=<ID={contig_id},length={contig_length}>\n")

                # Write column headers - apply name mapping for caller mode if provided
                if events:
                    sample_names = ["SAMPLE"]  # We'll use a single column for all samples
                else:
                    sample_names = ["SAMPLE"]

                header_line = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(sample_names) + "\n"
                f.write(header_line)

                for event in events:
                    # Prepare INFO field - apply name mapping if provided
                    info_items = []
                    for k, v in event.info.items():
                        if k == "SOURCES" and name_mapper:
                            # Convert file paths to display names
                            display_sources = name_mapper.convert_source_string(v)
                            info_items.append(f"SOURCES={display_sources}")
                        else:
                            info_items.append(f"{k}={v}")

                    info_field = ";".join(info_items)
                    if "SOURCES" not in info_field:
                        if name_mapper:
                            display_sources = name_mapper.convert_source_string(event.source_file)
                            info_field += f";SOURCES={display_sources}"
                        else:
                            info_field += f";SOURCES={event.source_file}"

                    # Get the FORMAT field from the representative event
                    format_field = event.format
                    format_keys = format_field.split(":")

                    # Prepare sample values
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

                        sample_part = "\t".join(sample_strings)
                    elif hasattr(event, "sample"):
                        formatted_values = self.format_sample_values(format_keys, event.sample)
                        if formatted_values.endswith(":.:."):
                            formatted_values = formatted_values[:-4]
                        sample_part = formatted_values
                    else:
                        sample_part = "./."

                    # Write the record
                    record_part1 = f"{event.chrom}\t{event.pos}\t{event.sv_id}\t{event.ref}\t{event.alt}\t"
                    record_part2 = f"{event.quality}\t{event.filter}\t{info_field}\t{format_field}\t"
                    f.write(record_part1 + record_part2 + sample_part + "\n")
