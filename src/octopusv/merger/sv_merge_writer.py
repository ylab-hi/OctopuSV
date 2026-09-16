import logging
import os

from octopusv.merger.name_mapper import default_label_from_path
from octopusv.utils.sample_consensus import resolve_sample_consensus
from octopusv.utils.source_path import normalize_source_path
from octopusv.utils.svcf_sample_parser import parse_svcf_sample_block
from octopusv.utils.vcf_info import format_vcf_info_item
from octopusv.utils.svcf_schema import (
    MODE_CALLER,
    mode_header,
    validate_source_id,
    validate_source_label,
    version_header,
)


class MergeWriterMixin:
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

        return ":".join(values)

    def _input_file_matches_event_sources(self, input_file, event_source_tokens):
        """Return whether an input file belongs to an event source set.

        event.source_file may contain full paths, basenames, or display-like
        stems. This helper accepts all of them.
        """
        input_file = str(input_file)
        input_basename = os.path.basename(input_file)
        input_stem = os.path.splitext(input_basename)[0]

        return (
            input_file in event_source_tokens
            or input_basename in event_source_tokens
            or input_stem in event_source_tokens
        )

    def _event_source_tokens(self, event):
        """Return normalized source tokens from event.source_file."""
        tokens = set()

        for raw_source in str(getattr(event, "source_file", "")).split(","):
            source = raw_source.strip()
            if not source:
                continue

            basename = os.path.basename(source)
            stem = os.path.splitext(basename)[0]

            tokens.add(source)
            tokens.add(basename)
            tokens.add(stem)

        return tokens

    def _sample_data_search_text(self, sample_name, sample_data):
        """Build searchable text from one sample/evidence block.

        In OctopuSV SVCF, the most useful source clues are often ID,
        original_id, and SC, for example:
            pbsv.INS.8
            svim.INS.7
            SVIM-v2.0.0
            Sniffles2.INS.5S0
        """
        parts = [str(sample_name)]

        if isinstance(sample_data, dict):
            for key in (
                "source_file",
                "original_source_file",
                "caller",
                "source",
                "original_id",
                "ID",
                "SC",
                "TY",
                "CO",
            ):
                value = sample_data.get(key)
                if value not in (None, "", "."):
                    parts.append(str(value))

            # Include the whole dict as a robust final fallback.
            parts.append(str(sample_data))
        else:
            parts.append(str(sample_data))

        return " ".join(parts).lower()

    def _infer_input_basename_for_sample_data(
        self,
        sample_name,
        sample_data,
        candidate_input_files,
        assigned_basenames,
    ):
        """Infer which input file one sample/evidence block came from.

        Priority:
            1. explicit source_file-like field
            2. evidence text containing input basename/stem
            3. sample_name containing input basename/stem
            4. no decision; caller will use a last-resort fallback
        """
        candidates = []
        for input_file in candidate_input_files:
            input_basename = os.path.basename(str(input_file))
            if input_basename not in assigned_basenames:
                candidates.append(str(input_file))

        if not candidates:
            return None

        # Method 1: explicit source_file field.
        if isinstance(sample_data, dict):
            source_file_field = str(sample_data.get("source_file", "")).lower()
            if source_file_field:
                for input_file in candidates:
                    input_basename = os.path.basename(input_file)
                    input_stem = os.path.splitext(input_basename)[0]

                    if (
                        input_file.lower() in source_file_field
                        or input_basename.lower() in source_file_field
                        or input_stem.lower() in source_file_field
                    ):
                        return input_basename

        search_text = self._sample_data_search_text(sample_name, sample_data)

        # Method 2: evidence content.
        # This fixes cases where merged_samples are internally ordered as
        # pbsv,svim even though the header order is sniffles,svim,pbsv.
        for input_file in candidates:
            input_basename = os.path.basename(input_file)
            input_stem = os.path.splitext(input_basename)[0].lower()

            if input_stem and input_stem in search_text:
                return input_basename

            if input_basename.lower() in search_text:
                return input_basename

        # Method 3: sample_name only. This is weaker than sample_data, but keeps
        # backward compatibility for older parsed objects.
        sample_name_text = str(sample_name).lower()
        for input_file in candidates:
            input_basename = os.path.basename(input_file)
            input_stem = os.path.splitext(input_basename)[0].lower()

            if input_stem and input_stem in sample_name_text:
                return input_basename

            if input_basename.lower() in sample_name_text:
                return input_basename

        return None

    @staticmethod
    def _normalize_source_path(source_file):
        """Return a stable key for exact input-file matching."""
        return normalize_source_path(source_file)

    def _display_name_for_input(self, input_file, name_mapper=None):
        """Return the output label for one input file."""
        if name_mapper is not None:
            return name_mapper.get_display_name(str(input_file))

        return default_label_from_path(input_file)

    @staticmethod
    def _source_id_from_sample_data(sample_data):
        """Return one original record ID from an evidence block."""
        if not isinstance(sample_data, dict):
            return "."

        # The evidence-block ID is the source record ID represented by this
        # block. ``original_id`` is retained only as a compatibility fallback
        # for older/manually-created objects that do not carry ID.
        source_id = sample_data.get(
            "ID",
            sample_data.get("original_id", "."),
        )

        if source_id in (None, "", "unknown"):
            return "."

        return str(source_id)

    def _prepare_caller_records_exact(self, event, name_mapper=None):
        """Build caller-mode records from exact source-to-evidence bindings.

        Fresh merge events carry ``merged_sample_records`` entries shaped as:
            (source_file, sample_name, sample_format, sample_data)

        Each entry remains a separate output evidence block. Therefore a source
        may appear more than once in SOURCES when it contributes more than one
        record to the same merged event. SOURCES, SOURCE_IDS, and evidence
        blocks are all generated from this same ordered list.

        Returns:
            list[dict] when exact source mapping is available.
            None when the event is a legacy object without merged_sample_records.

        Raises:
            ValueError if merged_sample_records exists but cannot be mapped
            unambiguously to the merge inputs.
        """
        records = getattr(event, "merged_sample_records", None)
        if records is None:
            return None

        if not all(
            isinstance(record, (tuple, list)) and len(record) == 4
            for record in records
        ):
            raise ValueError(
                "Invalid merged_sample_records: expected "
                "(source_file, sample_name, sample_format, sample_data) entries"
            )

        input_files = [str(path) for path in self.all_input_files]
        input_indices_by_path = {}

        for index, input_file in enumerate(input_files):
            normalized = self._normalize_source_path(input_file)
            input_indices_by_path.setdefault(normalized, []).append(index)

        ordered_records = []

        for record_order, record in enumerate(records):
            source_file, sample_name, sample_format, sample_data = record
            normalized_source = self._normalize_source_path(source_file)
            candidate_indices = input_indices_by_path.get(
                normalized_source,
                [],
            )

            if len(candidate_indices) != 1:
                if not candidate_indices:
                    reason = "does not match any merge input"
                else:
                    reason = "matches more than one merge input"

                raise ValueError(
                    "Cannot map caller evidence to exactly one input file: "
                    f"{source_file!r} {reason}."
                )

            input_index = candidate_indices[0]
            input_file = input_files[input_index]

            ordered_records.append(
                {
                    "input_index": input_index,
                    "record_order": record_order,
                    "source_name": self._display_name_for_input(
                        input_file,
                        name_mapper,
                    ),
                    "source_id": self._source_id_from_sample_data(
                        sample_data
                    ),
                    "sample_name": sample_name,
                    "sample_format": sample_format,
                    "sample_data": sample_data,
                }
            )

        ordered_records.sort(
            key=lambda record: (
                record["input_index"],
                record["record_order"],
            )
        )

        return ordered_records

    def _prepare_caller_records_legacy(self, event, name_mapper=None):
        """Best-effort compatibility path for old manually-created events.

        New merge output must use merged_sample_records instead. This fallback
        keeps older external objects usable, but source identity cannot always
        be determined exactly when that information was never stored.
        """
        self._legacy_caller_event_count = getattr(
            self,
            "_legacy_caller_event_count",
            0,
        ) + 1

        input_files = [str(path) for path in self.all_input_files]
        event_source_tokens = self._event_source_tokens(event)
        candidate_input_files = [
            input_file
            for input_file in input_files
            if self._input_file_matches_event_sources(
                input_file,
                event_source_tokens,
            )
        ]

        assigned_basenames = set()
        legacy_records = []

        for record_order, (
            sample_name,
            sample_format,
            sample_data,
        ) in enumerate(getattr(event, "merged_samples", [])):
            inferred_basename = self._infer_input_basename_for_sample_data(
                sample_name=sample_name,
                sample_data=sample_data,
                candidate_input_files=candidate_input_files,
                assigned_basenames=assigned_basenames,
            )

            if inferred_basename is None:
                inferred_basename = next(
                    (
                        os.path.basename(input_file)
                        for input_file in candidate_input_files
                        if os.path.basename(input_file) not in assigned_basenames
                    ),
                    None,
                )

            if inferred_basename is None:
                self._legacy_caller_unresolved_evidence_count = getattr(
                    self,
                    "_legacy_caller_unresolved_evidence_count",
                    0,
                ) + 1
                continue

            input_index = next(
                (
                    index
                    for index, input_file in enumerate(input_files)
                    if os.path.basename(input_file) == inferred_basename
                ),
                None,
            )

            if input_index is None:
                self._legacy_caller_unresolved_evidence_count = getattr(
                    self,
                    "_legacy_caller_unresolved_evidence_count",
                    0,
                ) + 1
                continue

            assigned_basenames.add(inferred_basename)
            input_file = input_files[input_index]

            legacy_records.append(
                {
                    "input_index": input_index,
                    "record_order": record_order,
                    "source_name": self._display_name_for_input(
                        input_file,
                        name_mapper,
                    ),
                    "source_id": self._source_id_from_sample_data(
                        sample_data
                    ),
                    "sample_name": sample_name,
                    "sample_format": sample_format,
                    "sample_data": sample_data,
                }
            )

        legacy_records.sort(
            key=lambda record: (
                record["input_index"],
                record["record_order"],
            )
        )

        return legacy_records

    def _prepare_caller_records(self, event, name_mapper=None):
        """Return caller evidence in deterministic input-file order."""
        exact_records = self._prepare_caller_records_exact(
            event,
            name_mapper=name_mapper,
        )

        if exact_records is not None:
            return exact_records

        return self._prepare_caller_records_legacy(
            event,
            name_mapper=name_mapper,
        )

    @staticmethod
    def _payload_from_sample_data(sample_data):
        if not isinstance(sample_data, dict):
            return None
        payload = sample_data.get("_octopusv_evidence_payload")
        return payload if isinstance(payload, dict) else None

    @staticmethod
    def _split_explicit_sources(source_value):
        """Split positional SOURCES without dropping placeholders."""
        if source_value in (None, ""):
            return []
        return [token.strip() for token in str(source_value).split(",")]

    def _caller_genotypes_from_payload(self, payload):
        """Return exact ``(caller, GT)`` pairs from one caller-mode payload.

        Source identity comes only from explicit SVCF data: INFO/SOURCES when
        present, otherwise the evidence block's SC field.  We never derive a
        caller from record ID, filename, or evidence order.
        """
        format_field = str(payload.get("format", "") or "")
        blocks = list(payload.get("blocks", ()) or ())
        if not format_field or not blocks:
            raise ValueError(
                "Sample consensus requires preserved caller evidence blocks."
            )

        parsed_blocks = [
            parse_svcf_sample_block(format_field, block)
            for block in blocks
        ]

        sources = self._split_explicit_sources(payload.get("sources"))
        if sources and len(sources) != len(parsed_blocks):
            raise ValueError(
                "Cannot synthesize sample consensus: SOURCES lists "
                f"{len(sources)} item(s), but the input record contains "
                f"{len(parsed_blocks)} caller evidence block(s)."
            )

        caller_genotypes = []
        for index, parsed in enumerate(parsed_blocks):
            source = sources[index] if sources else None
            if source in (None, "", "."):
                source = parsed.get("SC")

            if source in (None, "", ".", "OctopuSV"):
                raise ValueError(
                    "Cannot synthesize sample consensus because one caller "
                    "evidence block has no explicit caller identity."
                )

            caller_genotypes.append((str(source), parsed.get("GT", ".")))

        return caller_genotypes

    @staticmethod
    def _sample_length_from_record(record):
        svlen = record.get("svlen", ".")
        if svlen not in (None, "", "."):
            try:
                return str(abs(int(svlen)))
            except (TypeError, ValueError):
                pass
        return "."

    @staticmethod
    def _sample_coordinate_from_record(record):
        chrom = record.get("chrom", ".")
        pos = record.get("pos", ".")
        end_chrom = record.get("end_chrom", ".")
        end_pos = record.get("end_pos", ".")
        if any(value in (None, "", ".") for value in (chrom, pos, end_chrom, end_pos)):
            return "."
        return f"{chrom}_{pos}-{end_chrom}_{end_pos}"

    def _synthesize_sample_data(self, sample_records):
        """Synthesize one sample-mode column from exact caller evidence.

        ``sample_records`` contains one or more source records from the same
        biological-sample input that core merge placed in the same merged event.
        All underlying caller evidence contributes to genotype consensus, while
        record-level representation fields come from the first deterministic
        input record already chosen by the core merge ordering.  This function
        runs only after core grouping/representative selection, so it cannot
        change merge topology or representative-event selection.
        """
        if not sample_records:
            return None

        caller_genotypes = []
        payloads = []
        for sample_data in sample_records:
            payload = self._payload_from_sample_data(sample_data)
            if payload is None:
                raise ValueError(
                    "Sample-mode consensus requires preserved caller evidence "
                    "for every fresh merge record."
                )
            payloads.append(payload)
            caller_genotypes.extend(
                self._caller_genotypes_from_payload(payload)
            )

        consensus = resolve_sample_consensus(caller_genotypes)
        record = payloads[0].get("record")
        if not isinstance(record, dict):
            raise ValueError(
                "Sample-mode consensus is missing preserved record-level data."
            )

        return {
            "GT": consensus.gt,
            "AD": ".,.",
            "UC": str(consensus.unique_carrier_callers),
            "UV": str(consensus.valid_caller_votes),
            "LN": self._sample_length_from_record(record),
            "ST": str(record.get("strand", ".") or "."),
            "QV": str(record.get("quality", ".") or "."),
            "TY": str(record.get("svtype", ".") or "."),
            "ID": str(record.get("id", ".") or "."),
            "SC": "OctopuSV",
            "REF": str(record.get("ref", ".") or "."),
            "ALT": str(record.get("alt", ".") or "."),
            "CO": self._sample_coordinate_from_record(record),
        }

    def _prepare_events_for_sample_mode(self, events, name_mapper):
        """Prepare exact synthesized sample columns after core merge.

        Fresh merge events carry ``merged_sample_records`` entries shaped as::

            (source_file, sample_name, sample_format, sample_data)

        Exact source paths determine biological-sample columns.  For each input
        sample, every preserved caller evidence block from every contributing
        input record is reduced by the shared sample-consensus resolver.  No
        caller is selected merely because it appeared first.

        Legacy/manual events without exact bindings retain the historical
        best-effort mapping path and are reported in the returned summary.
        """
        processed_events = []
        input_files = [str(input_file) for input_file in self.all_input_files]

        input_indices_by_path = {}
        for index, input_file in enumerate(input_files):
            normalized_path = self._normalize_source_path(input_file)
            input_indices_by_path.setdefault(normalized_path, []).append(index)

        stats = {
            "legacy_sample_events": 0,
            "legacy_sample_unresolved_evidence_blocks": 0,
        }

        for event in events:
            merged_samples = list(getattr(event, "merged_samples", []))
            merged_sample_records = getattr(event, "merged_sample_records", None)

            source_to_records = {}
            unresolved_samples = []

            records_are_valid = (
                merged_sample_records is not None
                and len(merged_sample_records) == len(merged_samples)
                and all(
                    isinstance(record, (tuple, list)) and len(record) == 4
                    for record in merged_sample_records
                )
            )

            if records_are_valid:
                for record in merged_sample_records:
                    source_file, sample_name, sample_format, sample_data = record
                    normalized_source = self._normalize_source_path(source_file)
                    candidate_indices = input_indices_by_path.get(
                        normalized_source,
                        [],
                    )

                    if len(candidate_indices) != 1:
                        if not candidate_indices:
                            reason = "does not match any merge input"
                        else:
                            reason = "matches more than one merge input"
                        raise ValueError(
                            "Cannot map sample evidence to exactly one input file: "
                            f"{source_file!r} {reason}."
                        )

                    target_index = candidate_indices[0]
                    source_to_records.setdefault(target_index, []).append(sample_data)
            else:
                stats["legacy_sample_events"] += 1
                unresolved_samples.extend(merged_samples)

            source_to_sample = {
                target_index: self._synthesize_sample_data(sample_records)
                for target_index, sample_records in source_to_records.items()
            }

            if unresolved_samples:
                event_source_tokens = self._event_source_tokens(event)

                candidate_indices = [
                    index
                    for index, input_file in enumerate(input_files)
                    if self._input_file_matches_event_sources(
                        input_file,
                        event_source_tokens,
                    )
                ]
                candidate_input_files = [
                    input_files[index]
                    for index in candidate_indices
                ]

                for sample_name, sample_format, sample_data in unresolved_samples:
                    assigned_basenames = {
                        os.path.basename(input_files[index])
                        for index in source_to_sample
                    }

                    inferred_basename = self._infer_input_basename_for_sample_data(
                        sample_name=sample_name,
                        sample_data=sample_data,
                        candidate_input_files=candidate_input_files,
                        assigned_basenames=assigned_basenames,
                    )

                    target_index = None
                    if inferred_basename is not None:
                        for index in candidate_indices:
                            if index in source_to_sample:
                                continue
                            if os.path.basename(input_files[index]) == inferred_basename:
                                target_index = index
                                break

                    if target_index is None:
                        target_index = next(
                            (
                                index
                                for index in candidate_indices
                                if index not in source_to_sample
                            ),
                            None,
                        )

                    if target_index is not None:
                        source_to_sample[target_index] = sample_data
                    else:
                        stats[
                            "legacy_sample_unresolved_evidence_blocks"
                        ] += 1

            event.ordered_samples = [
                source_to_sample.get(index)
                for index in range(len(input_files))
            ]
            processed_events.append(event)

        return processed_events, stats

    def write_results(self, output_file, events, contigs, mode="caller", name_mapper=None, input_files=None):
        """Write merged results to output file with proper VCF header definitions.

        Args:
            output_file: Output file path
            events: List of merged events to write
            contigs: Dictionary of contig information
            mode: Merge mode ("caller" or "sample")
            name_mapper: NameMapper instance for name handling
            input_files: List of input file paths for dynamic header extraction
        """
        self._legacy_caller_event_count = 0
        self._legacy_caller_unresolved_evidence_count = 0

        if name_mapper and mode == "sample":
            processed_events, summary = self._prepare_events_for_sample_mode(
                events,
                name_mapper,
            )

            from .multi_sample_writer import MultiSampleWriter

            writer = MultiSampleWriter(name_mapper)
            writer.write_results(output_file, processed_events, contigs)

            legacy_events = summary.get("legacy_sample_events", 0)
            legacy_unresolved = summary.get(
                "legacy_sample_unresolved_evidence_blocks",
                0,
            )
            if legacy_events:
                logging.warning(
                    "Sample-mode writer used legacy source inference for %d "
                    "event(s); %d evidence block(s) could not be mapped.",
                    legacy_events,
                    legacy_unresolved,
                )

            return summary

        else:
            # Caller mode. SOURCES, SOURCE_IDS, and evidence blocks are all
            # generated from the same ordered source-to-record mapping.
            with open(output_file, "w") as f:
                self._write_vcf_header(
                    f,
                    contigs,
                    input_files,
                )

                for event in events:
                    ordered_records = self._prepare_caller_records(
                        event,
                        name_mapper=name_mapper,
                    )

                    sources_in_order = [
                        record["source_name"]
                        for record in ordered_records
                    ]
                    source_ids_in_order = [
                        record["source_id"]
                        for record in ordered_records
                    ]

                    # SVCF 1.1 positional INFO lists deliberately have no
                    # escaping layer.  Reject source labels/IDs that would
                    # make SOURCES or SOURCE_IDS ambiguous or syntactically
                    # invalid rather than silently emitting a corrupt file.
                    sources_in_order = [
                        validate_source_label(value)
                        for value in sources_in_order
                    ]
                    source_ids_in_order = [
                        validate_source_id(value)
                        for value in source_ids_in_order
                    ]

                    display_sources = (
                        ",".join(sources_in_order)
                        if sources_in_order
                        else "."
                    )
                    display_source_ids = (
                        ",".join(source_ids_in_order)
                        if source_ids_in_order
                        else "."
                    )

                    # Do not inherit stale merged source fields from the
                    # representative record. Rebuild them from this output's
                    # ordered evidence records.
                    info_items = []
                    for key, value in event.info.items():
                        if key in {"SOURCES", "SOURCE_IDS"}:
                            continue
                        info_items.append(format_vcf_info_item(key, value))

                    info_items.append(
                        f"SOURCES={display_sources}"
                    )
                    info_items.append(
                        f"SOURCE_IDS={display_source_ids}"
                    )
                    info_field = ";".join(info_items)

                    format_field = event.format
                    format_keys = format_field.split(":")

                    if ordered_records:
                        sample_strings = []
                        for record in ordered_records:
                            sample_data = record["sample_data"]

                            if isinstance(sample_data, dict):
                                sample_str = self.format_sample_values(
                                    format_keys,
                                    sample_data,
                                )
                            else:
                                sample_str = str(sample_data)

                            sample_strings.append(sample_str)

                        sample_part = "\t".join(sample_strings)
                    elif hasattr(event, "sample"):
                        sample_part = self.format_sample_values(
                            format_keys,
                            event.sample,
                        )
                    else:
                        sample_part = "./."

                    record_part1 = (
                        f"{event.chrom}\t{event.pos}\t{event.sv_id}\t"
                        f"{event.ref}\t{event.alt}\t"
                    )
                    record_part2 = (
                        f"{event.quality}\t{event.filter}\t{info_field}\t"
                        f"{format_field}\t"
                    )
                    f.write(
                        record_part1
                        + record_part2
                        + sample_part
                        + "\n"
                    )

            if self._legacy_caller_event_count:
                logging.warning(
                    "Caller-mode writer used legacy source inference for %d "
                    "event(s); %d evidence block(s) could not be mapped and "
                    "were omitted.",
                    self._legacy_caller_event_count,
                    self._legacy_caller_unresolved_evidence_count,
                )

            return {
                "legacy_caller_events": self._legacy_caller_event_count,
                "legacy_caller_unresolved_evidence_blocks": (
                    self._legacy_caller_unresolved_evidence_count
                ),
            }

    def _write_vcf_header(self, file_handle, contigs, input_files):
        """Write the caller-mode SVCF 1.1 header.

        Header generation is part of the SVCF contract. If dynamic header
        extraction fails, abort rather than appending a second fallback header
        to a partially written file.
        """
        # Extract and merge header definitions from all input files.
        merged_definitions = self._extract_and_merge_all_headers(input_files)

        import datetime

        file_handle.write("##fileformat=VCFv4.2\n")
        file_handle.write(version_header() + "\n")
        file_handle.write(mode_header(MODE_CALLER) + "\n")
        file_date = datetime.datetime.now().strftime("%Y-%m-%d|%I:%M:%S%p|")
        file_handle.write(f"##fileDate={file_date}\n")
        file_handle.write("##source=OctopuSV\n")

        for contig_id, contig_length in contigs.items():
            file_handle.write(
                f"##contig=<ID={contig_id},length={contig_length}>\n"
            )

        for alt_line in merged_definitions["alt_lines"]:
            file_handle.write(alt_line + "\n")

        # SOURCES / SOURCE_IDS are defined exactly once below.
        for info_line in merged_definitions["info_lines"]:
            info_id = self._extract_id_from_line(info_line, "INFO")
            if info_id in {"SOURCES", "SOURCE_IDS"}:
                continue
            file_handle.write(info_line + "\n")

        for filter_line in merged_definitions["filter_lines"]:
            file_handle.write(filter_line + "\n")

        for format_line in merged_definitions["format_lines"]:
            file_handle.write(format_line + "\n")

        file_handle.write(
            '##INFO=<ID=SOURCES,Number=.,Type=String,Description="Input source for each merged evidence block, in output order">\n'
        )
        file_handle.write(
            '##INFO=<ID=SOURCE_IDS,Number=.,Type=String,Description="Original record ID for each merged evidence block, in output order">\n'
        )

        file_handle.write(
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n"
        )

    def _extract_and_merge_all_headers(self, input_files):
        """Extract header definitions from all input files and merge them.

        Args:
            input_files: List of input file paths

        Returns:
            dict: Merged header definitions from all files
        """
        from octopusv.utils.svcf_utils import get_octopus_default_definitions, extract_original_header_definitions, \
            merge_header_definitions

        # Get OctopuSV default definitions
        octopus_defaults = get_octopus_default_definitions()

        # Collect all original definitions from all input files
        all_original_definitions = {
            'filter_lines': [],
            'info_lines': [],
            'format_lines': [],
            'alt_lines': []
        }

        # Track seen IDs to avoid duplicates
        seen_filter_ids = set()
        seen_info_ids = set()
        seen_format_ids = set()
        seen_alt_ids = set()

        # Extract definitions from each input file. Header extraction is part
        # of the output contract: if one input header cannot be read, abort
        # instead of silently dropping its custom definitions.
        if input_files:
            for input_file in input_files:
                file_headers = extract_original_header_definitions(input_file)

                # Add FILTER definitions (avoid duplicates)
                for line in file_headers.get('filter_lines', []):
                    filter_id = self._extract_id_from_line(line, 'FILTER')
                    if filter_id and filter_id not in seen_filter_ids:
                        all_original_definitions['filter_lines'].append(line)
                        seen_filter_ids.add(filter_id)

                # Add INFO definitions (avoid duplicates)
                for line in file_headers.get('info_lines', []):
                    info_id = self._extract_id_from_line(line, 'INFO')
                    if info_id and info_id not in seen_info_ids:
                        all_original_definitions['info_lines'].append(line)
                        seen_info_ids.add(info_id)

                # Add FORMAT definitions (avoid duplicates)
                for line in file_headers.get('format_lines', []):
                    format_id = self._extract_id_from_line(line, 'FORMAT')
                    if format_id and format_id not in seen_format_ids:
                        all_original_definitions['format_lines'].append(line)
                        seen_format_ids.add(format_id)

                # Add ALT definitions (avoid duplicates)
                for line in file_headers.get('alt_lines', []):
                    alt_id = self._extract_id_from_line(line, 'ALT')
                    if alt_id and alt_id not in seen_alt_ids:
                        all_original_definitions['alt_lines'].append(line)
                        seen_alt_ids.add(alt_id)

        # Merge with OctopuSV defaults
        return merge_header_definitions(all_original_definitions, octopus_defaults)

    def _extract_id_from_line(self, line, field_type):
        """Extract ID from header line.

        Args:
            line: Header line (e.g., ##INFO=<ID=SVTYPE,...)
            field_type: Type of field (INFO, FORMAT, FILTER, ALT)

        Returns:
            str: Extracted ID or None
        """
        try:
            if f'##{field_type}=<ID=' in line:
                return line.split('ID=')[1].split(',')[0]
        except (IndexError, AttributeError):
            pass
        return None

    def _write_basic_vcf_header(self, file_handle, contigs):
        """Write basic VCF header as fallback when dynamic generation fails.

        Args:
            file_handle: File handle to write to
            contigs: Dictionary of contig information
        """
        import datetime

        # Basic header information
        file_handle.write("##fileformat=VCFv4.2\n")
        file_handle.write(version_header() + "\n")
        file_handle.write(mode_header(MODE_CALLER) + "\n")
        file_date = datetime.datetime.now().strftime("%Y-%m-%d|%I:%M:%S%p|")
        file_handle.write(f"##fileDate={file_date}\n")
        file_handle.write("##source=OctopuSV\n")

        # Write contig information
        for contig_id, contig_length in contigs.items():
            file_handle.write(f"##contig=<ID={contig_id},length={contig_length}>\n")

        # Basic FILTER definitions
        file_handle.write('##FILTER=<ID=PASS,Description="All filters passed">\n')

        # Basic INFO definitions
        file_handle.write('##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">\n')
        file_handle.write(
            '##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant described in this record">\n')
        file_handle.write(
            '##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Difference in length between REF and ALT alleles">\n')
        file_handle.write('##INFO=<ID=CHR2,Number=1,Type=String,Description="Chromosome for end coordinate">\n')
        file_handle.write(
            '##INFO=<ID=SUPPORT,Number=1,Type=Integer,Description="Number of reads supporting this variant">\n')
        file_handle.write('##INFO=<ID=SVMETHOD,Number=1,Type=String,Description="Method used to detect SV">\n')
        file_handle.write('##INFO=<ID=STRAND,Number=1,Type=String,Description="Strand orientation of the SV">\n')
        file_handle.write(
            '##INFO=<ID=SOURCES,Number=.,Type=String,Description="Input source for each merged evidence block, in output order">\n')
        file_handle.write(
            '##INFO=<ID=SOURCE_IDS,Number=.,Type=String,Description="Original record ID for each merged evidence block, in output order">\n')

        # Basic ALT definitions
        file_handle.write('##ALT=<ID=DEL,Description="Deletion">\n')
        file_handle.write('##ALT=<ID=INV,Description="Inversion">\n')
        file_handle.write('##ALT=<ID=DUP,Description="Duplication">\n')
        file_handle.write('##ALT=<ID=INS,Description="Insertion">\n')
        file_handle.write('##ALT=<ID=TRA,Description="Translocation">\n')
        file_handle.write('##ALT=<ID=BND,Description="Breakend">\n')

        # Basic FORMAT definitions
        file_handle.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
        file_handle.write(
            '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">\n')
        file_handle.write('##FORMAT=<ID=LN,Number=1,Type=Integer,Description="Length of SV">\n')
        file_handle.write('##FORMAT=<ID=ST,Number=1,Type=String,Description="Strand orientation of SV">\n')
        file_handle.write('##FORMAT=<ID=QV,Number=1,Type=Integer,Description="Quality value">\n')
        file_handle.write('##FORMAT=<ID=TY,Number=1,Type=String,Description="Type of SV">\n')

        # Write column headers
        sample_names = ["SAMPLE"]
        header_line = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(sample_names) + "\n"
        file_handle.write(header_line)
