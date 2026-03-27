import logging
import re

from .base import Converter, get_alt_chrom_pos, get_bnd_pattern


class BNDPairToDELConverter(Converter):
    """Converter for identifying and merging BND pairs into DEL or INS events.

    DEL pattern: t[p[ event references a larger position than its own position
    Example: T[chr1:1357397[ at position 886721 references 1357397 (1357397 > 886721)

    INS detection: If the inserted sequence length >= 0.7 * abs(svLen),
    the event is classified as INS instead of DEL. This handles callers like
    GRIDSS that represent insertions as BND pairs with long inserted sequences.
    """

    def find_and_convert_pairs(self, events):
        """Find BND pairs and convert them to DEL or INS events.

        Args:
            events: List of BND events to process

        Returns:
            tuple: (remaining_events, converted_events)
        """
        remaining_events = []
        converted_events = []
        processed_events = set()

        # Use nested loop to check all pairs
        for i, event1 in enumerate(events):
            if id(event1) in processed_events:
                continue

            for j, event2 in enumerate(events[i + 1:], i + 1):
                if id(event2) in processed_events:
                    continue

                # Check if this pair forms a DEL or INS
                converted_event = self._check_and_convert_del_pair(event1, event2)
                if converted_event:
                    converted_events.append(converted_event)
                    processed_events.add(id(event1))
                    processed_events.add(id(event2))
                    sv_type = converted_event.info.get("SVTYPE", "unknown")
                    logging.debug(f"Converted BND pair {event1.id}-{event2.id} to {sv_type} event")
                    break

        # Add unprocessed events to remaining
        for event in events:
            if id(event) not in processed_events:
                remaining_events.append(event)

        return remaining_events, converted_events

    def _check_and_convert_del_pair(self, event1, event2):
        """Check if two events form a DEL or INS pair and convert them.

        DEL criteria:
        - One event has t[p[ pattern, the other has ]p]t pattern
        - They reference each other's positions
        - The t[p[ event references a position larger than its own position

        INS override:
        - If the inserted sequence length >= 0.7 * sv_len, classify as INS
        """
        try:
            # Both must be BND events
            if not (event1.is_BND() and event2.is_BND()):
                return None

            # Must be on same chromosome
            if event1.chrom != event2.chrom:
                return None

            pattern1 = get_bnd_pattern(event1.alt)
            pattern2 = get_bnd_pattern(event2.alt)

            chrom_alt1, pos_alt1 = get_alt_chrom_pos(event1.alt)
            chrom_alt2, pos_alt2 = get_alt_chrom_pos(event2.alt)

            # Check if they reference each other
            if not (event1.chrom == chrom_alt1 and event2.chrom == chrom_alt2 and
                    event1.pos == pos_alt2 and event2.pos == pos_alt1):
                return None

            # Check for DEL pattern: t[p[ references larger position
            if pattern1 == "t[p[" and pattern2 == "]p]t":
                # Check if t[p[ (event1) references a larger position
                if pos_alt1 > event1.pos:
                    start_pos = min(event1.pos, event2.pos)
                    end_pos = max(event1.pos, event2.pos)
                    base_event = event1 if event1.pos < event2.pos else event2
                    # Check if this is actually an INS (long inserted sequence, short distance)
                    ins_len = self._get_insertion_length(event1, event2)
                    sv_len = abs(end_pos - start_pos)
                    if sv_len > 0 and ins_len >= sv_len * 0.7:
                        return self._create_ins_event(base_event, start_pos, ins_len)
                    return self._create_del_event(base_event, start_pos, end_pos)

            elif pattern1 == "]p]t" and pattern2 == "t[p[":
                # Check if t[p[ (event2) references a larger position
                if pos_alt2 > event2.pos:
                    start_pos = min(event1.pos, event2.pos)
                    end_pos = max(event1.pos, event2.pos)
                    base_event = event1 if event1.pos < event2.pos else event2
                    # Check if this is actually an INS (long inserted sequence, short distance)
                    ins_len = self._get_insertion_length(event1, event2)
                    sv_len = abs(end_pos - start_pos)
                    if sv_len > 0 and ins_len >= sv_len * 0.7:
                        return self._create_ins_event(base_event, start_pos, ins_len)
                    return self._create_del_event(base_event, start_pos, end_pos)

            return None

        except Exception as e:
            logging.error(f"Error checking DEL/INS pair: {e}")
            return None

    def _get_insertion_length(self, event1, event2):
        """Extract insertion sequence length from ALT fields of a BND pair.

        The insertion sequence is the non-reference bases in the ALT field
        that are not part of the bracket notation.
        For example: GATTTTGTTC[chr1:876082[ has insertion 'ATTTTGTTC' (9bp).
        We take the maximum insertion length from both events in the pair.
        """
        max_ins = 0
        for event in [event1, event2]:
            alt = event.alt
            # Remove bracket and coordinate part, e.g. [chr1:12345[ or ]chr1:12345]
            seq = re.sub(r'[\[\]][^\[\]]+:\d+[\[\]]', '', alt)
            # The remaining is ref base + inserted sequence
            # Subtract 1 for the ref base
            ins_len = max(0, len(seq) - 1)
            if ins_len > max_ins:
                max_ins = ins_len
        return max_ins

    def _create_del_event(self, base_event, start_pos, end_pos):
        """Create a DEL event from base BND event."""
        import copy
        del_event = copy.deepcopy(base_event)

        # Set to start position
        del_event.pos = start_pos

        # Modify to DEL format
        del_event.alt = "<DEL>"
        del_event.info["SVTYPE"] = "DEL"
        del_event.info["END"] = end_pos
        del_event.info["SVLEN"] = abs(end_pos - start_pos)
        del_event.info["CHR2"] = del_event.chrom
        del_event.info["SVMETHOD"] = "OctopuSV"

        return del_event

    def _create_ins_event(self, base_event, start_pos, ins_len):
        """Create an INS event from base BND event.

        For INS events, END = POS + SVLEN following the SVCF convention
        used by OctopuSV for interval-based comparison.
        """
        import copy
        ins_event = copy.deepcopy(base_event)

        ins_event.pos = start_pos
        ins_event.alt = "<INS>"
        ins_event.info["SVTYPE"] = "INS"
        ins_event.info["END"] = str(start_pos + ins_len)
        ins_event.info["SVLEN"] = str(ins_len)
        ins_event.info["CHR2"] = ins_event.chrom
        ins_event.info["SVMETHOD"] = "OctopuSV"

        return ins_event

    def convert(self, event):
        """Standard convert method for compatibility with base Converter class.

        Note: This method is not used for pair processing.
        Use find_and_convert_pairs() instead for batch processing.
        """
        pass