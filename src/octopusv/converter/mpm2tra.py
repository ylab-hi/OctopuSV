from .base import (
    Converter,
    compare_chromosomes,
    get_alt_chrom_pos,
    is_same_bnd_event,
)


def merged_mate_support(first, second):
    """SUPPORT of one junction written as two mate records.

    Both ends describe the same junction and the same reads, so support is not
    added.  With per-sample values from one file (same samples, same order)
    the larger value is taken in each sample and then summed; otherwise the
    larger known record value is kept.  Unknown stays unknown.
    """
    a = getattr(first, "support_per_sample", None)
    b = getattr(second, "support_per_sample", None)
    if a is not None and b is not None and len(a) == len(b):
        if any(x is None and y is None for x, y in zip(a, b)):
            return "."
        return str(sum(max(v for v in (x, y) if v is not None) for x, y in zip(a, b)))
    known = []
    for event in (first, second):
        value = str(event.info.get("SUPPORT", "."))
        if value.isdigit():
            known.append(int(value))
    return str(max(known)) if known else "."


class MatePairMergeToTRAConverter(Converter):
    """Class to transform a pair of identical BND events into a single TRA event."""

    def convert(self, pair):
        event1, event2 = pair
        # Check if the pair of events has the same BND pattern
        if is_same_bnd_event(event1, event2):
            # Determine which event should be retained
            retained_event = event1 if compare_chromosomes(event1, event2) else event2
            other_event = event2 if retained_event is event1 else event1

            retained_event.info["SUPPORT"] = merged_mate_support(retained_event, other_event)

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
