from .base import (
    Converter,
    compare_chromosomes,
    get_alt_chrom_pos,
    is_same_bnd_event,
)


class MatePairMergeToTRAConverter(Converter):
    """Class to transform a pair of identical BND events into a single TRA event."""

    def convert(self, pair):
        event1, event2 = pair
        # Check if the pair of events has the same BND pattern
        if is_same_bnd_event(event1, event2):
            # Determine which event should be retained
            retained_event = event1 if compare_chromosomes(event1, event2) else event2

            # Modify the SVTYPE to TRA
            retained_event.info["SVTYPE"] = "TRA"
            retained_event.info["SVLEN"] = "."
            chrom_alt, pos_alt = get_alt_chrom_pos(retained_event.alt)
            retained_event.info["END"] = pos_alt
            retained_event.info["CHR2"] = chrom_alt
            retained_event.info["SVMETHOD"] = "OctopuSV"

            # Return the retained event only
            return [retained_event]

        # If the pair of events does not have the same BND pattern, return an empty list
        return []
