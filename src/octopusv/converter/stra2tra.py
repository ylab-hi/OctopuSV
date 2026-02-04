from .base import Converter, get_alt_chrom_pos


class SingleTRAToTRAConverter(Converter):
    """This class inherits from the `Converter` base class and implements the conversion logic for single TRA events."""

    def convert(self, event):
        event.info["SVTYPE"] = "TRA"
        event.info["SVLEN"] = "."
        chrom_alt, pos_alt = get_alt_chrom_pos(event.alt)
        event.info["END"] = pos_alt
        event.info["CHR2"] = chrom_alt
        event.info["SVMETHOD"] = "OctopuSV"
