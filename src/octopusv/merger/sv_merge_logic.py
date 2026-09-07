def should_merge(event1, event2, max_distance=None, max_length_ratio=None, min_jaccard=0.10):
    """Determine whether two ordinary SV events should be merged.

    By default, OctopuSV keeps the v0.4.2 adaptive distance and
    SV-type-specific length-ratio thresholds. Explicit max_distance or
    max_length_ratio values override those adaptive defaults.

    Interval Jaccard is applied only to DEL, DUP, and INV. INS continues
    to use breakpoint distance and length ratio only.

    Args:
        event1: First SV event.
        event2: Second SV event.
        max_distance: Optional global breakpoint-distance override. If None,
            use the SV-type- and size-aware adaptive threshold.
        max_length_ratio: Optional global SV-length-ratio override. If None,
            use the SV-type-specific threshold.
        min_jaccard: Minimum interval Jaccard for DEL, DUP, and INV.
            Set to 0 to disable the Jaccard requirement.

    Returns:
        bool: True if the events should be merged.
    """
    # Compare chromosomes.
    if event1.chrom != event2.chrom:
        return False

    # Get event type.
    sv_type = event1.info.get("SVTYPE", "")

    # Safely get lengths.
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

    # Use explicit user overrides when provided; otherwise preserve the
    # v0.4.2 adaptive/type-specific defaults.
    if max_distance is None:
        distance_threshold = _get_distance_threshold(
            sv_type, min(length1, length2)
        )
    else:
        distance_threshold = max_distance

    if max_length_ratio is None:
        length_ratio_threshold = _get_length_ratio_threshold(sv_type)
    else:
        length_ratio_threshold = max_length_ratio

    # Breakpoint comparison.
    start_diff = abs(event1.start_pos - event2.start_pos)
    end_diff = abs(event1.end_pos - event2.end_pos)

    if start_diff > distance_threshold:
        return False

    # INS is a point-like event on the reference and does not use END
    # compatibility or interval Jaccard.
    if sv_type != "INS" and end_diff > distance_threshold:
        return False

    # Length-ratio comparison.
    if length1 <= 0 or length2 <= 0:
        return False

    ratio = max(length1, length2) / min(length1, length2)
    if ratio > length_ratio_threshold:
        return False

    # DEL, DUP, and INV describe reference intervals, so require a minimum
    # interval overlap. Jaccard=0 disables this extra constraint and
    # reproduces the v0.4.2 ordinary-SV matching logic when no other
    # overrides are supplied.
    if sv_type in {"DEL", "DUP", "INV"} and min_jaccard > 0:
        jaccard = _calculate_interval_jaccard(event1, event2)
        if jaccard < min_jaccard:
            return False

    return True


def _calculate_interval_jaccard(event1, event2):
    """Calculate inclusive interval Jaccard for two SV events."""
    intersection = max(
        0,
        min(event1.end_pos, event2.end_pos)
        - max(event1.start_pos, event2.start_pos)
        + 1,
    )

    union = (
        max(event1.end_pos, event2.end_pos)
        - min(event1.start_pos, event2.start_pos)
        + 1
    )

    if union <= 0:
        return 0.0

    return intersection / union


def _get_distance_threshold(sv_type, min_length):
    """Get the v0.4.2 adaptive distance threshold."""
    base_threshold = {
        "INS": 200,  # More lenient for insertions
        "DEL": 150,
        "DUP": 100,
        "INV": 100,
    }.get(sv_type, 50)

    # Adjust threshold based on event size.
    if min_length >= 1000:
        return base_threshold * 2
    if min_length >= 500:
        return base_threshold * 1.5
    return base_threshold


def _get_length_ratio_threshold(sv_type):
    """Get the v0.4.2 SV-type-specific length-ratio threshold."""
    return {
        "INS": 3.0,  # Very lenient for insertions
        "DEL": 2.0,  # Somewhat lenient for deletions
        "DUP": 1.5,
        "INV": 1.5,
    }.get(sv_type, 1.3)
