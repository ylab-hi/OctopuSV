from .base import Converter, get_alt_chrom_pos, get_bnd_pattern


def is_mate_pair_reciprocal_translocation(event1, event2):
    # Define the set of qualified pairings of BND patterns, used by MatePairReciprocalTranslocationToTRAConverter
    qualified_pairings = [
        {"]p]t"},
        {"t[p["},
        {"t]p]", "[p[t"},
    ]
    # Extract the patterns from each event
    event1_pattern = get_bnd_pattern(event1.alt)
    event2_pattern = get_bnd_pattern(event2.alt)

    # Check if the patterns match one of the qualified pairings
    return any({event1_pattern, event2_pattern} == pairing for pairing in qualified_pairings)


class MatePairReciprocalTranslocationToTRAConverter(Converter):
    """This class inherits from the `Converter` base class and implements.

    the conversion logic for reciprocal translocation.
    """

    def convert(self, pair):
        """Check if this pair satisfies the criteria for reciprocal translocation."""
        event1, event2 = pair
        if is_mate_pair_reciprocal_translocation(event1, event2):
            # Convert events
            self.convert_to_TRA(event1, event2)
            return [event1, event2]  # Return a list of transformed events

        return []  # If the pair doesn't satisfy the criteria, return an empty list

    def convert_to_TRA(self, event1, event2):
        """Convert a pair of events to reciprocal translocation."""
        # Modify event1
        event1.info["SVTYPE"] = "TRA"
        chrom_alt, pos_alt = get_alt_chrom_pos(event1.alt)
        event1.info["END"] = pos_alt
        event1.info["CHR2"] = chrom_alt
        event1.info["RTID"] = event2.id
        event1.info["SVLEN"] = "."
        event1.info["SVMETHOD"] = "OctopuSV"
        # Modify event2
        event2.info["SVTYPE"] = "TRA"
        chrom_alt, pos_alt = get_alt_chrom_pos(event2.alt)
        event2.info["END"] = pos_alt
        event2.info["CHR2"] = chrom_alt
        event2.info["RTID"] = event1.id
        event2.info["SVLEN"] = "."
        event2.info["SVMETHOD"] = "OctopuSV"
        # No need to return anything as the states of event1 and event2 are modified in-place
