from .base import (
    Converter,
    get_alt_chrom_pos,
    is_SpecialNoMateDiffBndPairReciprocalTranslocation,
)


class SpecialNoMateDiffBNDPairReciprocalTranslocationToTRAConverter(Converter):
    """This class inherits from the `Converter` base class and implements the conversion logic for reciprocal translocation."""

    def convert(self, pair):
        event1, event2 = pair
        # Check if this pair satisfies the criteria for reciprocal translocation
        if is_SpecialNoMateDiffBndPairReciprocalTranslocation(event1, event2):
            # Convert events
            self.convert_to_TRA(event1, event2)
            return [event1, event2]  # Return a list of transformed events
        return []  # If the pair doesn't satisfy the criteria, return an empty list

    def convert_to_TRA(self, event1, event2):
        # Convert a pair of events to reciprocal translocation
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
