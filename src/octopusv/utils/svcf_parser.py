import logging
import os
import re


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
        format (str): Format of the sample data related to the SV event.
        sample (dict): Sample-specific data for the SV event.
        source_file (str): The source file from which this SV event was parsed.
        sv_type (str): The SV type of the event.
        bnd_pattern (str): The BND pattern from the ALT field, if applicable.
        start_chrom (str): Start chromosome of the SV.
        start_pos (int): Start position of the SV.
        end_chrom (str): End chromosome of the SV.
        end_pos (int): End position of the SV.
        sample_name (str): Name of the sample.
    """

    def __init__(self, chrom, pos, sv_id, ref, alt, quality, filter, info, format, sample, source_file, sample_name):
        self.chrom = chrom
        self.pos = int(pos)
        self.sv_id = sv_id
        self.ref = ref
        self.alt = alt
        self.quality = quality
        self.filter = filter
        self.info = self._parse_info(info)
        self.format = format
        self.sample = self._parse_sample(sample)
        self.source_file = source_file
        self.sample_name = sample_name
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
        """Parses the sample string based on the format to extract all relevant data.
        Handles fields that may contain colons (like ALT and CO).
        """
        format_keys = self.format.split(":")
        sample_parts = sample.split(":")
        result = {}

        # Special handling of fields that may contain colons
        special_fields = ["ALT", "CO", "REF"]  # Fields that might need special handling
        special_field_index = None

        for i, key in enumerate(format_keys):
            if key in special_fields:
                special_field_index = i
                break

        if special_field_index is not None:
            # Regular fields
            for i in range(special_field_index):
                result[format_keys[i]] = sample_parts[i]

            # Process the special field and all fields after it
            result[format_keys[special_field_index]] = ":".join(sample_parts[special_field_index:])
        else:
            # If there are no special fields, handle them in the usual way
            result = dict(zip(format_keys, sample_parts, strict=False))

        return result

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
        # For other SV types
        co = self.sample.get("CO")
        if co and "-" in co:
            start, end = co.split("-")
            start_chrom, start_pos = start.split("_")
            end_chrom, end_pos = end.split("_")
            return start_chrom, int(start_pos), end_chrom, int(end_pos)
        start_chrom = self.chrom
        start_pos = self.pos
        end_chrom = self.info.get("CHR2", self.chrom)
        end_pos_str = self.info.get("END", self.pos)
        if end_pos_str == "." or end_pos_str is None:
            # Try to use SVLEN
            svlen = self.info.get("SVLEN", None)
            if svlen and svlen != ".":
                try:
                    end_pos = self.pos + abs(int(svlen))
                except ValueError:
                    end_pos = self.pos
            else:
                end_pos = self.pos
        else:
            try:
                end_pos = int(end_pos_str)
            except ValueError:
                logging.warning(
                    f"Invalid end_pos '{end_pos_str}' for SV {self.sv_id} in {self.source_file}, setting end_pos to start_pos."
                )
                end_pos = self.pos
        return start_chrom, start_pos, end_chrom, end_pos


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
        """Parses each file in the filenames list, creating SVCFEvent objects for each valid SV record."""
        for filename in self.filenames:
            with open(filename) as file:
                sample_name = None
                for line in file:
                    if line.startswith("#CHROM"):
                        header_parts = line.strip().split("\t")
                        if len(header_parts) > 9:
                            sample_name = header_parts[9]
                        else:
                            sample_name = os.path.basename(filename)  # Use filename as sample name if not provided
                        continue
                    if line.startswith("#"):
                        continue  # Skip comment lines that start with '#'
                    parts = line.strip().split("\t")
                    if len(parts) < 10:
                        continue  # Ensure that there are enough parts to form a complete SV event
                    sv_event = SVCFEvent(
                        *parts[:10], source_file=filename, sample_name=sample_name
                    )  # Create an SV event object
                    self.events.append(sv_event)  # Add the event to the list of events
