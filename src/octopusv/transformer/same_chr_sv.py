from .base import EventTransformer
from octopusv.converter.bnd2del import BNDPairToDELConverter
from octopusv.converter.bnd2dup_pair import BNDPairToDUPConverter
from octopusv.converter.bnd2inv_pair import BNDPairToINVConverter


class SameChrSVTransformer(EventTransformer):
    """Class for transforming SV events on the same chromosome.

    This transformer handles:
    1. BND pair merging (DEL, DUP, INV)
    2. Remaining BND events are kept as BND (not forced to TRA)
    """

    def apply_transforms(self, events):
        """Apply all transformation strategies to a list of events."""

        # First, try to convert DEL pairs
        del_converter = BNDPairToDELConverter()
        events_after_del, converted_del_events = del_converter.find_and_convert_pairs(events)

        # Then, try to convert DUP pairs from remaining events
        dup_converter = BNDPairToDUPConverter()
        events_after_dup, converted_dup_events = dup_converter.find_and_convert_pairs(events_after_del)

        # Then, try to convert INV pairs from remaining events
        inv_converter = BNDPairToINVConverter()
        events_to_process, converted_inv_events = inv_converter.find_and_convert_pairs(events_after_dup)

        # Finally, apply regular converters to remaining events
        for event in events_to_process:
            for strategy in self.transform_strategies:
                strategy.convert(event)

        # Combine all events
        return events_to_process + converted_del_events + converted_dup_events + converted_inv_events

    def write_vcf(self, headers, events, output_file):
        """Write the transformed events to a VCF file."""
        with open(output_file, "w") as f:
            for header in headers:
                f.write(header + "\n")
            for event in events:
                f.write(str(event) + "\n")