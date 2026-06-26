from .bnd_merge_logic import parse_bnd_alt


_PATTERN_TABLE = {
    "t[p[": ("left", "right", False),
    "t]p]": ("left", "left", True),
    "]p]t": ("right", "left", False),
    "[p[t": ("right", "right", True),
}


def should_merge_tra(event1, event2, delta=100, min_overlap_ratio=0.4, strand_consistency=True):
    """Determine whether two TRA events should be merged.

    TRA/BND records are only safe to merge when they describe the same
    breakpoint cluster AND, when bracket ALT is available for both records, the
    same physical breakend adjacency.

    This avoids over-merging caller evidence such as:
        A: t[B[   (A-left joined to B-right, non-RC)
        A: [B[t   (A-right joined to B-right, RC)
    which are coordinate-near but topology-incompatible.

    Reciprocal mate representations are still allowed:
        A: t[B[   ==   B: ]A]t
    """
    tra_delta = delta * 2

    if {event1.start_chrom, event1.end_chrom} != {event2.start_chrom, event2.end_chrom}:
        return False

    if strand_consistency:
        sig1 = _breakend_orientation_signature(event1)
        sig2 = _breakend_orientation_signature(event2)
        if sig1 is not None and sig2 is not None:
            return _same_breakend_adjacency(sig1, sig2, tra_delta)

    return _legacy_coordinate_match(event1, event2, tra_delta)


def _breakend_orientation_signature(event):
    """Return a side-labeled physical adjacency signature for a TRA/BND event."""
    pattern, target_chr, target_pos = parse_bnd_alt(getattr(event, "alt", None))
    if pattern not in _PATTERN_TABLE:
        return None

    local_side, remote_side, reverse_complement = _PATTERN_TABLE[pattern]

    local_chrom = getattr(event, "start_chrom", None) or getattr(event, "chrom", None)
    local_pos = getattr(event, "start_pos", None) or getattr(event, "pos", None)

    remote_chrom = target_chr or getattr(event, "end_chrom", None)
    remote_pos = target_pos or getattr(event, "end_pos", None)

    if local_chrom is None or local_pos is None or remote_chrom is None or remote_pos is None:
        return None

    return (
        (str(local_chrom), int(local_pos), local_side),
        (str(remote_chrom), int(remote_pos), remote_side),
        bool(reverse_complement),
    )


def _same_endpoint(a, b, tol):
    return (
        a[0] == b[0]
        and a[2] == b[2]
        and abs(int(a[1]) - int(b[1])) <= tol
    )


def _same_breakend_adjacency(sig1, sig2, tol):
    """Compare physical BND adjacencies after reciprocal-mate normalization."""
    a_local, a_remote, a_rc = sig1
    b_local, b_remote, b_rc = sig2

    if a_rc != b_rc:
        return False

    same_frame = (
        _same_endpoint(a_local, b_local, tol)
        and _same_endpoint(a_remote, b_remote, tol)
    )

    reciprocal_frame = (
        _same_endpoint(a_local, b_remote, tol)
        and _same_endpoint(a_remote, b_local, tol)
    )

    return same_frame or reciprocal_frame


def _legacy_coordinate_match(event1, event2, tra_delta):
    if _should_swap_positions(event1.alt, event2.alt):
        e2_start = event2.end_pos
        e2_end = event2.start_pos
    else:
        e2_start = event2.start_pos
        e2_end = event2.end_pos

    start_diff = abs(event1.start_pos - e2_start)
    end_diff = abs(event1.end_pos - e2_end)

    return not (start_diff > tra_delta or end_diff > tra_delta)


def _should_swap_positions(alt1, alt2):
    """Determine if positions should be swapped based on BND patterns.

    This is retained only for the legacy coordinate fallback path.
    """
    if not alt1 or not alt2:
        return False

    pattern1 = _classify_bnd_pattern(alt1)
    pattern2 = _classify_bnd_pattern(alt2)

    return _are_reciprocal_patterns(pattern1, pattern2)


def _classify_bnd_pattern(alt):
    """Classify BND pattern from ALT field for the legacy coordinate fallback."""
    if not alt:
        return "UNKNOWN"

    import re

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
    reciprocal_pairs = {("RIGHT_TO_LEFT", "LEFT_TO_RIGHT"), ("RIGHT_TO_RIGHT", "LEFT_TO_LEFT")}
    return (pattern1, pattern2) in reciprocal_pairs or (pattern2, pattern1) in reciprocal_pairs
