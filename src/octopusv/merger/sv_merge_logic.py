def should_merge(
    event1,
    event2,
    max_distance=None,
    max_length_ratio=None,
    min_jaccard=0.10,
):
    """Determine whether two ordinary SV events should be merged.

    If max_distance or max_length_ratio is not explicitly provided,
    OctopuSV preserves the v0.4.2 adaptive/type-specific defaults.

    Interval Jaccard is applied only to DEL, DUP, and INV.
    INS does not use interval Jaccard.
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

    # Use explicit user overrides when provided.
    # Otherwise preserve the v0.4.2 defaults.
    if max_distance is None:
        distance_threshold = _get_distance_threshold(
            sv_type,
            min(length1, length2),
        )
    else:
        distance_threshold = max_distance

    if max_length_ratio is None:
        length_ratio_threshold = _get_length_ratio_threshold(sv_type)
    else:
        length_ratio_threshold = max_length_ratio

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
    if ratio > length_ratio_threshold:
        return False

    # DEL/DUP/INV represent reference intervals.
    # Apply Jaccard only to these SV types.
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
    """Get distance threshold based on SV type and minimum length."""
    base_threshold = {
        "INS": 200,
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


def get_max_distance_threshold(sv_type):
    """Return the maximum adaptive distance threshold for an SV type.

    This public helper is used by the optimized active-window pruning in
    sv_merger.py. It must remain at least as permissive as should_merge()
    when no explicit --max-distance override is supplied.
    """
    return _get_distance_threshold(sv_type, 1000)


def _get_length_ratio_threshold(sv_type):
    """Get length ratio threshold based on SV type."""
    return {
        "INS": 3.0,
        "DEL": 2.0,
        "DUP": 1.5,
        "INV": 1.5,
    }.get(sv_type, 1.3)
