"""Observed member breakpoint offsets in a finished OctopuSV event group.

These are spreads of *reported member coordinates*, NOT caller confidence
intervals, physical breakpoint uncertainty, or additional merge criteria.
No change to grouping, representative selection or SVCF coordinate semantics.
"""

from __future__ import annotations

MEMBER_POS_RANGE = "MERGE_POS_RANGE"
MEMBER_END_RANGE = "MERGE_END_RANGE"
MEMBER_RANGE_KEYS = frozenset((MEMBER_POS_RANGE, MEMBER_END_RANGE))
MEMBER_RANGE_HEADERS = {
    MEMBER_POS_RANGE: (
        '##INFO=<ID=MERGE_POS_RANGE,Number=2,Type=Integer,Description="Minimum and maximum offsets '
        'of merged input member POS coordinates relative to representative POS; '
        'reported coordinates, not a confidence interval. '
        'For TRA/BND, member endpoints are aligned to representative chromosome sides. '
        'In sample mode members are input sample-level SV records.">'
    ),
    MEMBER_END_RANGE: (
        '##INFO=<ID=MERGE_END_RANGE,Number=2,Type=Integer,Description="Minimum and maximum offsets '
        'of merged input member second-breakpoint coordinates relative to representative '
        'second breakpoint; reported coordinates, not a confidence interval. '
        'For TRA/BND, member endpoints are aligned to representative chromosome sides. '
        'Not emitted for insertions (INS END is length-derived). '
        'In sample mode members are input sample-level SV records.">'
    ),
}


def needs_member_range_headers(events):
    """Range keys present in events, in stable header order."""
    seen = set()
    for event in events:
        seen.update(MEMBER_RANGE_KEYS.intersection(getattr(event, "info", {})))
    return tuple(key for key in (MEMBER_POS_RANGE, MEMBER_END_RANGE) if key in seen)


def _coordinates(event):
    """Physical endpoints from parsed SVCFEvent, with chromosome identity."""
    try:
        return (
            str(event.start_chrom), int(event.start_pos),
            str(event.end_chrom), int(event.end_pos),
        )
    except (AttributeError, TypeError, ValueError):
        return None


def _aligned_coordinates(member, representative):
    """Match both endpoints to the same chromosome sides as representative."""
    m = _coordinates(member)
    r = _coordinates(representative)
    if m is None or r is None:
        return None
    chrom_a, pos_a, chrom_b, pos_b = m
    ref_a, _, ref_b, _ = r
    if (chrom_a, chrom_b) == (ref_a, ref_b):
        return pos_a, pos_b
    # Only interchromosomal TRA endpoints can be unambiguously flipped using
    # chromosome identity. Same-chromosome endpoints are left in input order.
    if ref_a != ref_b and (chrom_b, chrom_a) == (ref_a, ref_b):
        return pos_b, pos_a
    return None


def _range_text(offsets):
    return f"{min(offsets)},{max(offsets)}"


def add_member_breakpoint_ranges(representative, members):
    """Replace inherited range fields with this group's actual member offsets.

    In sample mode the inputs are already sample-level records: this measures
    differences among their representative breakpoints, not all original
    caller positions. Singletons deliberately have no range fields.
    """
    for key in MEMBER_RANGE_KEYS:
        representative.info.pop(key, None)
    if len(members) < 2:
        return
    svtype = str(representative.info.get("SVTYPE", getattr(representative, "sv_type", "")))
    if svtype not in {"DEL", "DUP", "INV", "INS", "TRA", "BND"}:
        return
    ref = _coordinates(representative)
    if ref is None:
        return

    offsets_pos = []
    offsets_end = []
    for member in members:
        aligned = _aligned_coordinates(member, representative)
        if aligned is None:
            # One unresolved or incompatible member: do not publish a partial
            # range falsely suggesting it summarizes every member.
            return
        offsets_pos.append(aligned[0] - ref[1])
        if svtype != "INS":
            offsets_end.append(aligned[1] - ref[3])

    representative.info[MEMBER_POS_RANGE] = _range_text(offsets_pos)
    if svtype != "INS":
        representative.info[MEMBER_END_RANGE] = _range_text(offsets_end)
