import logging
import os
import re

from .svcf_sample_parser import parse_svcf_sample_block
from .svcf_coordinate_parser import parse_svcf_co
from .text_io import open_text_auto


class SVCFEvent:
    """Represents a structural variant (SV) event parsed from an SVCF file.

    Attributes:
        chrom (str): Chromosome on which the SV event occurs.
        pos (int): Starting position of the SV on the chromosome.
        sv_id (str): Unique identifier for the SV event.
        ref (str): Reference bases at the SV position.
        alt (str): Alternate bases indicating the SV.
        quality (str): Quality score of the SV event.
        filter (str): Filter status of the SV event.
        info (dict): Dictionary containing additional information about the SV event.
        raw_info (str): The ORIGINAL, unparsed INFO string. Kept verbatim so that
            inspect/validate can detect duplicate INFO keys (info dict collapses
            repeated keys to the last value).
        format (str): Format of the sample data related to the SV event.
        sample (dict): Sample-specific data for the SV event (first sample).
        raw_sample_columns (list): Raw string of EVERY sample column (multi-sample safe).
        sample_names (list): All trailing column names from the #CHROM header,
            one per sample column. Sample-mode needs the full list to map each
            column to its sample identity; falls back to [sample_name].
        source_file (str): The source file from which this SV event was parsed.
        sv_type (str): The SV type of the event.
        bnd_pattern (str): The BND pattern from the ALT field, if applicable.
        start_chrom (str): Start chromosome of the SV.
        start_pos (int): Start position of the SV.
        end_chrom (str): End chromosome of the SV.
        end_pos (int): End position of the SV.
        sample_name (str): Name of the sample.
    """

    # Matches a CO coordinate token like "chr1_10002-chr1_10056".
    _CO_PATTERN = re.compile(r"([^\s:]+_\d+-[^\s:]+_\d+)")

    def __init__(self, chrom, pos, sv_id, ref, alt, quality, filter, info, format, sample,
                 source_file, sample_name, raw_sample_columns=None, sample_names=None):
        self.chrom = chrom
        self.pos = int(pos)
        self.sv_id = sv_id
        self.ref = ref
        self.alt = alt
        self.quality = quality
        self.filter = filter
        # Keep the raw INFO string BEFORE parsing into a dict. _parse_info()
        # collapses duplicate keys (e.g. a buggy merge writing SOURCE_IDS twice),
        # so the only way for inspect/validate to see duplicates is the raw text.
        self.raw_info = info
        self.info = self._parse_info(info)
        self.format = format
        self.sample = self._parse_sample(sample)

        # Multi-sample support: keep EVERY raw sample column string, one per sample.
        # Single-sample input -> [sample]; multi-sample -> [col1, col2, ...].
        # Mirrors SVEvent.samples so both event classes share one mental model.
        # NOTE: assigned BEFORE _parse_coordinates(), which now consults it.
        self.raw_sample_columns = list(raw_sample_columns) if raw_sample_columns is not None else [sample]

        self.source_file = source_file
        self.sample_name = sample_name
        # All trailing sample/column names from the #CHROM header. Sample-mode
        # inspect aligns each raw_sample_columns entry to its sample identity by
        # this list. Falls back to [sample_name] when only one column exists.
        self.sample_names = list(sample_names) if sample_names is not None else [sample_name]
        self.sv_type = self.info.get("SVTYPE", "")
        self.bnd_pattern = self._extract_bnd_pattern()

        try:
            self.start_chrom, self.start_pos, self.end_chrom, self.end_pos = self._parse_coordinates()
        except Exception as e:
            logging.warning(f"Warning: Error parsing coordinates for {self.sv_id} in {self.source_file}: {e}")
            # Use default values as fallback options
            self.start_chrom, self.start_pos = self.chrom, self.pos
            self.end_chrom = self.info.get("CHR2", self.chrom)
            end_value = self.info.get("END", self.pos)
            if end_value == "." or end_value is None:
                self.end_pos = self.start_pos
            else:
                try:
                    self.end_pos = int(end_value)
                except ValueError:
                    logging.error(f"Invalid END value: {end_value}, setting end_pos to start_pos")
                    self.end_pos = self.start_pos

    def _extract_bnd_pattern(self):
        """Extract the breakend pattern from ALT field if present."""
        if self.sv_type in ("BND", "TRA"):
            return self.alt
        return None

    def _parse_info(self, info_str):
        """Parses the info field from a VCF record into a dictionary.

        Args:
            info_str (str): The raw info string from a VCF record.

        Returns:
            dict: A dictionary where each key is an info field name and the value is the field value.
        """
        info = {}
        for item in info_str.split(";"):
            if "=" in item:
                key, value = item.split("=", 1)  # Only split at the first '='
                info[key] = value
            else:
                info[item] = True  # Flags without a value are stored as True
        return info

    def _parse_sample(self, sample):
        """Parse one SVCF evidence block using the shared SVCF parser."""
        result = parse_svcf_sample_block(self.format, sample)

        # Keep the record-level ID used by existing merge/sample-mode code.
        # The evidence-level ID remains available as result["ID"].
        result["original_id"] = self.sv_id

        return result

    def _find_co_with_coords(self):
        """Return the first structurally valid CO coordinate token.

        The shared parser is conservative and understands contig names that
        contain underscores or hyphens. Ambiguous CO text is never guessed.
        """
        for col in self.raw_sample_columns:
            co_val = col.rsplit(":", 1)[-1]
            if parse_svcf_co(co_val) is not None:
                return co_val

        fallback = self.sample.get("CO")
        return fallback if parse_svcf_co(fallback) is not None else None

    def _parse_coordinates(self):
        """Extracts and parses the coordinates of the SV from the ALT field or INFO field."""
        if self.sv_type in ("BND", "TRA"):
            alt = self.alt
            # Define regex pattern to match ALT field for BND/TRA events
            pattern = re.compile(r"([ACGTNacgtn]*)([\[\]])([^:\[\]]+):(\d+)([\[\]])([ACGTNacgtn]*)")
            match = pattern.match(alt)
            if match:
                seq_before, bracket1, end_chrom, end_pos_str, bracket2, seq_after = match.groups()
                try:
                    end_pos = int(end_pos_str)
                except ValueError:
                    logging.warning(
                        f"Invalid end_pos '{end_pos_str}' for SV {self.sv_id} in {self.source_file}, setting end_pos to start_pos."
                    )
                    end_pos = self.pos
                return self.chrom, self.pos, end_chrom, end_pos
            # If ALT field parsing fails, try to get coordinates from INFO
            end_chrom = self.info.get("CHR2", self.chrom)
            end_pos_str = self.info.get("END", self.pos)
            try:
                end_pos = int(end_pos_str)
            except ValueError:
                logging.warning(
                    f"Invalid end_pos '{end_pos_str}' for SV {self.sv_id} in {self.source_file}, setting end_pos to start_pos."
                )
                end_pos = self.pos
            return self.chrom, self.pos, end_chrom, end_pos

        # For ordinary SVs, the record-level CHROM/POS/INFO-END fields are
        # authoritative. FORMAT/CO belongs to per-sample/per-caller evidence and
        # must not override the coordinates of a merged representative record.
        start_chrom = self.chrom
        start_pos = self.pos
        end_chrom = self.info.get("CHR2", self.chrom)
        end_pos_str = self.info.get("END")

        if end_pos_str not in (None, "", "."):
            try:
                return start_chrom, start_pos, end_chrom, int(end_pos_str)
            except (ValueError, TypeError):
                logging.warning(
                    f"Invalid end_pos '{end_pos_str}' for SV {self.sv_id} in {self.source_file}; "
                    "trying record-level SVLEN before CO fallback."
                )

        # If END is absent or unusable, prefer record-level SVLEN over FORMAT/CO.
        # This keeps a merged representative tied to its own POS/SVLEN rather than
        # borrowing coordinates from another caller's evidence block.
        svlen = self.info.get("SVLEN")
        if svlen not in (None, "", "."):
            try:
                end_pos = self.pos + abs(int(svlen))
                return start_chrom, start_pos, end_chrom, end_pos
            except (ValueError, TypeError):
                pass

        # Legacy fallback: recover an end coordinate from CO only when the
        # record-level END and SVLEN are both unavailable or unusable. Preserve
        # the record-level CHROM/POS even in this fallback.
        co = self._find_co_with_coords()
        parsed_co = self._coords_from_co(co)
        if parsed_co is not None:
            _co_start_chrom, _co_start_pos, co_end_chrom, co_end_pos = parsed_co
            return start_chrom, start_pos, co_end_chrom, co_end_pos

        return start_chrom, start_pos, end_chrom, start_pos

    def _coords_from_co(self, co):
        """Parse CO through the shared unambiguous coordinate parser."""
        return parse_svcf_co(co)



class SVCFFileEventCreator:
    """Parses all SVCF files to create SVCFEvent instances.

    Attributes:
        filenames (list): List of filenames (strings) to be parsed.
        events (list): List of SVCFEvent objects parsed from the files.
    """

    def __init__(self, filenames):
        self.filenames = filenames
        self.events = []

    def parse(self):
        """Parse all SVCF records, failing loudly on malformed data rows."""
        for filename in self.filenames:
            with open_text_auto(filename) as file:
                sample_name = None
                sample_names = None
                for line_number, line in enumerate(file, 1):
                    if line.startswith("#CHROM"):
                        header_parts = line.strip().split("\t")
                        if len(header_parts) > 9:
                            sample_names = header_parts[9:]
                            sample_name = sample_names[0]
                        else:
                            sample_name = os.path.basename(filename)
                            sample_names = [sample_name]
                        continue
                    if line.startswith("#") or not line.strip():
                        continue

                    parts = line.strip().split("\t")
                    if len(parts) < 10:
                        raise ValueError(
                            f"Malformed SVCF record in {filename!r} on line "
                            f"{line_number}: expected at least 10 tab-separated "
                            f"columns, got {len(parts)}."
                        )

                    if sample_name is None:
                        raise ValueError(
                            f"Malformed SVCF {filename!r}: data record on line "
                            f"{line_number} appears before a #CHROM header."
                        )

                    sv_event = SVCFEvent(
                        *parts[:10], source_file=filename, sample_name=sample_name,
                        raw_sample_columns=parts[9:], sample_names=sample_names,
                    )
                    self.events.append(sv_event)
