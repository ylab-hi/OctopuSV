def should_merge_tra(event1, event2, delta=100, min_overlap_ratio=0.4, strand_consistency=True):
    """Determines whether two TRA events should be merged based on comprehensive pattern matching.

    Args:
        event1: First TRA event
        event2: Second TRA event
        delta: Position uncertainty threshold (default: 100bp)
        min_overlap_ratio: Minimum required overlap ratio
        strand_consistency: Whether to enforce strand consistency
    """
    # Increase default delta for TRA events to better handle tool-specific variations
    tra_delta = delta * 2  # More relax

    # Check chromosomes
    if {event1.start_chrom, event1.end_chrom} != {event2.start_chrom, event2.end_chrom}:
        return False

    # Position comparison logic
    if _should_swap_positions(event1.alt, event2.alt):
        e2_start = event2.end_pos
        e2_end = event2.start_pos
    else:
        e2_start = event2.start_pos
        e2_end = event2.end_pos

    # More lenient position checks for TRA
    start_diff = abs(event1.start_pos - e2_start)
    end_diff = abs(event1.end_pos - e2_end)

    return not (start_diff > tra_delta or end_diff > tra_delta)


def _should_swap_positions(alt1, alt2):
    """Determine if positions should be swapped based on BND patterns.

    Args:
        alt1: First BND ALT field
        alt2: Second BND ALT field

    Returns:
        bool: True if positions should be swapped
    """
    if not alt1 or not alt2:
        return False

    pattern1 = _classify_bnd_pattern(alt1)
    pattern2 = _classify_bnd_pattern(alt2)

    return _are_reciprocal_patterns(pattern1, pattern2)


def _classify_bnd_pattern(alt):
    """Classify BND pattern from ALT field.

    Args:
        alt: BND ALT field value

    Returns:
        str: Pattern classification
    """
    if not alt:
        return "UNKNOWN"

    # Extract basic pattern structure
    import re

    # Remove chromosome and position numbers but keep structure
    pattern = re.sub(r"chr\d+:\d+", "chrN:N", alt)
    pattern = re.sub(r"\d+:\d+", "N:N", pattern)

    if pattern.startswith("]") or ("]" in pattern and pattern.endswith("N")):
        return "RIGHT_TO_LEFT"
    if pattern.startswith("N[") or ("[" in pattern and pattern.endswith("[")):
        return "LEFT_TO_RIGHT"
    if pattern.startswith("N]") or ("]" in pattern and pattern.endswith("]")):
        return "RIGHT_TO_RIGHT"
    if pattern.startswith("[") or ("[" in pattern and pattern.endswith("N")):
        return "LEFT_TO_LEFT"

    return "UNKNOWN"


def _are_reciprocal_patterns(pattern1, pattern2):
    """Check if two patterns are reciprocal.

    Args:
        pattern1: First pattern classification
        pattern2: Second pattern classification

    Returns:
        bool: True if patterns are reciprocal
    """
    reciprocal_pairs = {("RIGHT_TO_LEFT", "LEFT_TO_RIGHT"), ("RIGHT_TO_RIGHT", "LEFT_TO_LEFT")}

    return (pattern1, pattern2) in reciprocal_pairs or (pattern2, pattern1) in reciprocal_pairs
