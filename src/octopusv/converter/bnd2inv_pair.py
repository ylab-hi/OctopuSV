import logging
import re
from .base import Converter, get_alt_chrom_pos, get_bnd_pattern


class BNDPairToINVConverter(Converter):
    """Converter for identifying and merging BND pairs into INV events.

    INV pattern: Both events have the same pattern (either both ]p]t or both [p[t)
    Example: A]chr1:10978252] and T]chr1:10978048] -> INV
    """

    def find_and_convert_pairs(self, events):
        """Find BND pairs and convert them to INV events.

        Args:
            events: List of BND events to process

        Returns:
            tuple: (remaining_events, converted_inv_events)
        """
        remaining_events = []
        converted_inv_events = []
        processed_events = set()

        # Use nested loop to check all pairs
        for i, event1 in enumerate(events):
            if id(event1) in processed_events:
                continue

            for j, event2 in enumerate(events[i + 1:], i + 1):
                if id(event2) in processed_events:
                    continue

                # Check if this pair forms an INV
                inv_event = self._check_and_convert_inv_pair(event1, event2)
                if inv_event:
                    converted_inv_events.append(inv_event)
                    processed_events.add(id(event1))
                    processed_events.add(id(event2))
                    logging.debug(f"Converted BND pair {event1.id}-{event2.id} to INV event")
                    break

        # Add unprocessed events to remaining
        for event in events:
            if id(event) not in processed_events:
                remaining_events.append(event)

        return remaining_events, converted_inv_events

    def _check_and_convert_inv_pair(self, event1, event2):
        """Check if two events form an INV pair and convert them if so.

        INV criteria:
        - Both events have the same BND pattern (]p]t or [p[t)
        - They reference each other's positions
        - Same chromosome
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

            # Check for INV pattern: both events have same pattern
            if ((pattern1 == "t]p]" and pattern2 == "t]p]") or
                    (pattern1 == "[p[t" and pattern2 == "[p[t")):
                # Create INV event from the one with smaller position
                start_pos = min(event1.pos, event2.pos)
                end_pos = max(event1.pos, event2.pos)
                base_event = event1 if event1.pos < event2.pos else event2
                return self._create_inv_event(base_event, start_pos, end_pos)

            return None

        except Exception as e:
            logging.error(f"Error checking INV pair: {e}")
            return None

    def _create_inv_event(self, base_event, start_pos, end_pos):
        """Create an INV event from base BND event."""
        import copy
        inv_event = copy.deepcopy(base_event)

        # Set to start position
        inv_event.pos = start_pos

        # Modify to INV format
        inv_event.alt = "<INV>"
        inv_event.info["SVTYPE"] = "INV"
        inv_event.info["END"] = end_pos
        inv_event.info["SVLEN"] = abs(end_pos - start_pos)
        inv_event.info["CHR2"] = inv_event.chrom
        inv_event.info["SVMETHOD"] = "OctopuSV"

        return inv_event

    def convert(self, event):
        """Standard convert method for compatibility with base Converter class.

        Note: This method is not used for pair processing.
        Use find_and_convert_pairs() instead for batch processing.
        """
        pass