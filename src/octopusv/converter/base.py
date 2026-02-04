import logging
import re

from natsort import natsorted

logging.basicConfig(level=logging.INFO)


class Converter:
    """This is an abstract base class for all converter classes.
    It provides a common interface for all converters.
    The `convert` method is a placeholder that needs to be overridden in each concrete converter class.
    """

    def convert(self, event):
        """This method should be overridden in a subclass."""
        raise NotImplementedError


def get_bnd_pattern(alt):
    """Get the pattern of BND: t[p[, t]p], ]p]t, [p[t.
    This function uses regular expressions to match the ALT field,
    considering possible multiple nucleotides in the sequences and special chromosome names.
    """
    patterns = {
        "t[p[": r"^[ACGTNacgtn]*\[[^\[\]]+:\d+\[$",
        "t]p]": r"^[ACGTNacgtn]*\][^\[\]]+:\d+\]$",
        "[p[t": r"^\[[^\[\]]+:\d+\[[ACGTNacgtn]*$",
        "]p]t": r"^\][^\[\]]+:\d+\][ACGTNacgtn]*$",
    }

    for pattern_name, regex in patterns.items():
        if re.match(regex, alt):
            return pattern_name
    return None


def get_alt_chrom_pos(alt):
    """Extract chromosome and position from ALT field in VCF file.
    This function handles both standard chromosome names and special ones like hs37d5.
    """
    try:
        # Handle various BND format patterns
        brackets = {"[", "]"}
        if not any(b in alt for b in brackets):
            return None, None

        parts = re.split(r"[\[\]]", alt)
        for part in parts:
            if ":" in part:
                chrom, pos = part.split(":")
                try:
                    return chrom, int(pos)
                except ValueError:
                    continue

        logging.info(f"Unexpected ALT format, it should be something like N]chr10:69650962] or with sequences: {alt}")
        return None, None

    except Exception as e:
        logging.info(f"Error parsing ALT field: {alt}, error: {e!s}")
        return None, None


def is_same_bnd_event(event1, event2) -> bool:
    """Define whether the BND represent the same TRA events, used by MatePairMergeToTRAConverter."""
    qualified_pairings = [
        {"t[p[", "t]p]"},
        {"]p]t", "t[p["},  # Added pattern.
    ]
    event1_pattern = get_bnd_pattern(event1.alt)
    event2_pattern = get_bnd_pattern(event2.alt)
    return any({event1_pattern, event2_pattern} == pairing for pairing in qualified_pairings)


def is_independent_bnd_event(event1, event2):
    """Determine if two mate_pair_bnd events are independent events."""
    qualified_pairings = [
        {"t]p]", "]p]t"},
        {"t]p]"},
        {"[p[t"},
        {"t[p[", "t]p]"},
        {"t[p[", "[p[t"},
        {"]p]t", "[p[t"},
    ]
    event1_pattern = get_bnd_pattern(event1.alt)
    event2_pattern = get_bnd_pattern(event2.alt)
    return any({event1_pattern, event2_pattern} == pairing for pairing in qualified_pairings)


def is_independent_special_bnd_event(event1, event2):
    """Determine if two special_no_mate_diff_bnd_pair events are independent events."""
    qualified_pairings = [
        {"t[p[", "t]p]"},
        {"t[p[", "[p[t"},
        {"t]p]", "]p]t"},
        {"]p]t", "[p[t"},
    ]
    event1_pattern = get_bnd_pattern(event1.alt)
    event2_pattern = get_bnd_pattern(event2.alt)
    return any({event1_pattern, event2_pattern} == pairing for pairing in qualified_pairings)


def is_SpecialNoMateDiffBndPairReciprocalTranslocation(event1, event2):
    """Determine if SpecialNoMateDiffBndPairReciprocalTranslocation."""
    qualified_pairings = [{"t[p[", "]p]t"}, {"t]p]", "[p[t"}]
    event1_pattern = get_bnd_pattern(event1.alt)
    event2_pattern = get_bnd_pattern(event2.alt)
    return any({event1_pattern, event2_pattern} == pairing for pairing in qualified_pairings)


def compare_chromosomes(event1, event2):
    """Determine if the chromosome of event1 is smaller than the chromosome of event2 using natural sort order.

    Returns:
        bool: True if the chromosome of event1 is smaller or equal to the chromosome of event2, False otherwise.
    """
    original_order = [event1.chrom, event2.chrom]
    # Get the order when the chromosomes are naturally sorted
    sorted_order = natsorted(original_order)
    # Return True if the original order matches the sorted order
    return original_order == sorted_order
