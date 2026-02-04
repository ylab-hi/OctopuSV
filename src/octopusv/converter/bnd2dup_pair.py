import logging
import re
from .base import Converter, get_alt_chrom_pos, get_bnd_pattern


class BNDPairToDUPConverter(Converter):
    """Converter for identifying and merging BND pairs into DUP events.

    DUP pattern: t[p[ event references a smaller position than its own position
    Example: C[chr1:10004[ at position 10574 references 10004 (10004 < 10574)
    """

    def find_and_convert_pairs(self, events):
        """Find BND pairs and convert them to DUP events.

        Args:
            events: List of BND events to process

        Returns:
            tuple: (remaining_events, converted_dup_events)
        """
        remaining_events = []
        converted_dup_events = []
        processed_events = set()

        # Use nested loop to check all pairs
        for i, event1 in enumerate(events):
            if id(event1) in processed_events:
                continue

            for j, event2 in enumerate(events[i + 1:], i + 1):
                if id(event2) in processed_events:
                    continue

                # Check if this pair forms a DUP
                dup_event = self._check_and_convert_dup_pair(event1, event2)
                if dup_event:
                    converted_dup_events.append(dup_event)
                    processed_events.add(id(event1))
                    processed_events.add(id(event2))
                    logging.debug(f"Converted BND pair {event1.id}-{event2.id} to DUP event")
                    break

        # Add unprocessed events to remaining
        for event in events:
            if id(event) not in processed_events:
                remaining_events.append(event)

        return remaining_events, converted_dup_events

    def _check_and_convert_dup_pair(self, event1, event2):
        """Check if two events form a DUP pair and convert them if so.

        DUP criteria:
        - One event has t[p[ pattern, the other has ]p]t pattern
        - They reference each other's positions
        - The t[p[ event references a position smaller than its own position
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

            # Check for DUP pattern: t[p[ references smaller position
            if pattern1 == "t[p[" and pattern2 == "]p]t":
                # Check if t[p[ (event1) references a smaller position
                if pos_alt1 < event1.pos:
                    # This is a DUP pattern
                    start_pos = min(event1.pos, event2.pos)
                    end_pos = max(event1.pos, event2.pos)
                    # Use the event at smaller position as base
                    base_event = event1 if event1.pos < event2.pos else event2
                    return self._create_dup_event(base_event, start_pos, end_pos)

            elif pattern1 == "]p]t" and pattern2 == "t[p[":
                # Check if t[p[ (event2) references a smaller position
                if pos_alt2 < event2.pos:
                    # This is a DUP pattern
                    start_pos = min(event1.pos, event2.pos)
                    end_pos = max(event1.pos, event2.pos)
                    # Use the event at smaller position as base
                    base_event = event1 if event1.pos < event2.pos else event2
                    return self._create_dup_event(base_event, start_pos, end_pos)

            return None

        except Exception as e:
            logging.error(f"Error checking DUP pair: {e}")
            return None

    def _create_dup_event(self, base_event, start_pos, end_pos):
        """Create a DUP event from base BND event."""
        import copy
        dup_event = copy.deepcopy(base_event)

        # Set to start position
        dup_event.pos = start_pos

        # Modify to DUP format
        dup_event.alt = "<DUP>"
        dup_event.info["SVTYPE"] = "DUP"
        dup_event.info["END"] = end_pos
        dup_event.info["SVLEN"] = abs(end_pos - start_pos)
        dup_event.info["CHR2"] = dup_event.chrom
        dup_event.info["SVMETHOD"] = "OctopuSV"

        return dup_event

    def convert(self, event):
        """Standard convert method for compatibility with base Converter class.

        Note: This method is not used for pair processing.
        Use find_and_convert_pairs() instead for batch processing.
        """
        pass