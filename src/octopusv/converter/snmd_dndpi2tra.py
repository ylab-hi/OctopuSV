from .base import Converter, get_alt_chrom_pos, is_independent_special_bnd_event


class SpecialNoMateDiffBNDPairIndependentToTRAConverter(Converter):
    """This class implements the conversion logic for independent translocation."""

    def convert(self, pair):
        event1, event2 = pair
        # Check if this pair satisfies the criteria for independent translocation
        if is_independent_special_bnd_event(event1, event2):
            # Convert events
            self.convert_to_TRA(event1, event2)
            return [event1, event2]  # Return a list of transformed events

        return []  # If the pair doesn't satisfy the criteria, return an empty list

    def convert_to_TRA(self, event1, event2):
        # Convert a pair of independent events to TRA
        # Modify event1
        event1.info["SVTYPE"] = "TRA"
        event1.info["SVLEN"] = "."
        chrom_alt, pos_alt = get_alt_chrom_pos(event1.alt)
        event1.info["END"] = pos_alt
        event1.info["CHR2"] = chrom_alt
        event1.info["SVMETHOD"] = "OctopuSV"
        # Modify event2
        event2.info["SVTYPE"] = "TRA"
        event2.info["SVLEN"] = "."
        chrom_alt, pos_alt = get_alt_chrom_pos(event2.alt)
        event2.info["END"] = pos_alt
        event2.info["CHR2"] = chrom_alt
        event2.info["SVMETHOD"] = "OctopuSV"
