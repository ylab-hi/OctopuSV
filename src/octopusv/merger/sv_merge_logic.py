def should_merge(event1, event2, max_distance=50, max_length_ratio=1.3, min_jaccard=0.7):
    """Determines whether two SV events should be merged based on simple position and length comparison.

    Args:
        event1: First SV event
        event2: Second SV event
        max_distance: Base maximum allowed distance between positions
        max_length_ratio: Base maximum allowed ratio between lengths
        min_jaccard: Base minimum required Jaccard index (not used in simplified version)

    Returns:
        bool: True if events should be merged
    """
    # Compare chromosomes
    if event1.chrom != event2.chrom:
        return False

    # Get event type
    sv_type = event1.info.get("SVTYPE", "")

    # Safely get lengths
    def get_length(event):
        svlen = event.info.get("SVLEN", ".")
        try:
            if svlen == "." or not svlen:
                return event.end_pos - event.start_pos + 1
            return abs(int(svlen))
        except (ValueError, TypeError):
            return event.end_pos - event.start_pos + 1

    length1 = get_length(event1)
    length2 = get_length(event2)

    # Get type-specific thresholds
    distance_threshold = _get_distance_threshold(sv_type, min(length1, length2))
    length_ratio_threshold = _get_length_ratio_threshold(sv_type)

    # Position comparison
    start_diff = abs(event1.start_pos - event2.start_pos)
    end_diff = abs(event1.end_pos - event2.end_pos)

    if start_diff > distance_threshold:
        return False

    # For non-INS events, also check end position
    if sv_type != "INS" and end_diff > distance_threshold:
        return False

    # Length ratio comparison
    if length1 == 0 or length2 == 0:
        return False

    ratio = max(length1, length2) / min(length1, length2)
    return ratio <= length_ratio_threshold


def _get_distance_threshold(sv_type, min_length):
    """Get distance threshold based on SV type and minimum length."""
    base_threshold = {
        "INS": 200,  # More lenient for insertions
        "DEL": 150,
        "DUP": 100,
        "INV": 100,
    }.get(sv_type, 50)

    # Adjust threshold based on event size
    if min_length >= 1000:
        return base_threshold * 2
    if min_length >= 500:
        return base_threshold * 1.5
    return base_threshold


def _get_length_ratio_threshold(sv_type):
    """Get length ratio threshold based on SV type."""
    return {
        "INS": 3.0,  # Very lenient for insertions
        "DEL": 2.0,  # Somewhat lenient for deletions
        "DUP": 1.5,
        "INV": 1.5,
    }.get(sv_type, 1.3)
