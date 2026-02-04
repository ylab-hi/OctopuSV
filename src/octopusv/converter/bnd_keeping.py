import logging
from .base import Converter, get_alt_chrom_pos, get_bnd_pattern


class BNDKeepingConverter(Converter):
    """This class inherits from the Converter base class and keeps BND events as BND.

    For remaining BND events that were not paired or converted to other SV types,
    we keep them as BND rather than forcing conversion to TRA, as suggested by reviewers.
    These events may represent genuine single breakpoints, complex rearrangements,
    or incomplete detections that should remain classified as BND.
    """

    def convert(self, event):
        """Keep BND events as BND with updated info fields."""
        try:
            if event.is_BND():
                # Simply update the info fields but keep SVTYPE as BND
                self.update_bnd_info(event)
        except Exception as e:
            logging.error("Failed to update BND event info: ", e)

    def update_bnd_info(self, event):
        """Update BND event info fields while preserving BND type."""
        chrom_alt, pos_alt = get_alt_chrom_pos(event.alt)

        # Update basic info fields
        if chrom_alt and pos_alt:
            event.info["END"] = pos_alt
            event.info["CHR2"] = chrom_alt
        else:
            event.info["END"] = "."
            event.info["CHR2"] = "."

        # Keep SVTYPE as BND (don't change it)
        event.info["SVTYPE"] = "BND"
        event.info["SVLEN"] = "."  # BND events don't have meaningful SVLEN
        event.info["SVMETHOD"] = "OctopuSV"