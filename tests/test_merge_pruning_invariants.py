from __future__ import annotations

from types import SimpleNamespace

import pytest

from octopusv.merger.bnd_merger import BNDMerger
from octopusv.merger.sv_merger import SVMerger
from octopusv.merger.tra_merger import TRAMerger


def _group_signature(groups):
    """Compare membership only; representative choice/order is not science here."""
    return tuple(
        sorted(
            tuple(sorted(event.sv_id for event in group))
            for group in groups
        )
    )


def _ordinary_event(svtype: str, event_id: str, start: int, end: int):
    return SimpleNamespace(
        chrom="chr1",
        start_pos=start,
        end_pos=end,
        sv_id=event_id,
        info={"SVTYPE": svtype, "SVLEN": str(abs(end - start + 1))},
    )


def _ordinary_specs(svtype: str):
    if svtype == "INS":
        # INS ignores end-coordinate distance in should_merge().  Include both
        # the short-event default threshold and the long-event adaptive maximum.
        return [
            ("short.seed", 100, 149),
            ("short.near", 110, 164),
            ("short.boundary", 300, 349),       # +200: merge boundary
            ("short.outside", 301, 350),        # +201: outside short threshold
            ("long.seed", 10_000, 10_999),
            ("long.boundary", 10_400, 11_399),  # +400: adaptive max boundary
            ("long.outside", 10_401, 11_400),   # +401: outside adaptive max
        ]

    base = {"DEL": 150, "DUP": 100, "INV": 100}[svtype]
    adaptive_max = {"DEL": 300, "DUP": 200, "INV": 200}[svtype]
    return [
        ("short.seed", 100, 199),
        ("short.overlap", 110, 209),
        ("short.boundary", 100 + base, 199 + base),
        ("short.outside", 101 + base, 200 + base),
        ("long.seed", 10_000, 10_999),
        ("long.boundary", 10_000 + adaptive_max, 10_999 + adaptive_max),
        ("long.outside", 10_001 + adaptive_max, 11_000 + adaptive_max),
    ]


def _optimized_ordinary_groups(svtype: str, min_jaccard: float):
    events = [_ordinary_event(svtype, *spec) for spec in _ordinary_specs(svtype)]
    classified = {svtype: {"chr1": events}}
    merger = SVMerger(
        classified,
        all_input_files=[],
        min_jaccard=min_jaccard,
    )
    merger.merge()
    return merger.event_groups[svtype]["chr1"]


def _full_scan_ordinary_groups(svtype: str, min_jaccard: float):
    """Use production first-match grouping, but call its full-scan entry point.

    This intentionally does not reimplement should_merge().  The only difference
    from merge() is that add_and_merge_event() scans every existing group rather
    than using the active-window pruning path.
    """
    events = [_ordinary_event(svtype, *spec) for spec in _ordinary_specs(svtype)]
    merger = SVMerger(
        classified_events={},
        all_input_files=[],
        min_jaccard=min_jaccard,
    )
    merger.merged_events = {svtype: {"chr1": []}}
    merger.event_groups = {svtype: {"chr1": []}}

    for event in sorted(events, key=lambda e: (e.start_pos, e.end_pos, e.sv_id)):
        merger.add_and_merge_event(svtype, "chr1", event)

    return merger.event_groups[svtype]["chr1"]


@pytest.mark.parametrize("svtype", ["DEL", "DUP", "INV", "INS"])
@pytest.mark.parametrize("min_jaccard", [0.0, 0.05, 0.10])
def test_ordinary_active_window_preserves_full_scan_group_membership(
    svtype,
    min_jaccard,
):
    optimized = _group_signature(
        _optimized_ordinary_groups(svtype, min_jaccard)
    )
    full_scan = _group_signature(
        _full_scan_ordinary_groups(svtype, min_jaccard)
    )

    assert optimized == full_scan


def _tra_event(event_id, start_chrom, start_pos, end_chrom, end_pos, alt):
    return SimpleNamespace(
        sv_id=event_id,
        chrom=start_chrom,
        pos=start_pos,
        start_chrom=start_chrom,
        start_pos=start_pos,
        end_chrom=end_chrom,
        end_pos=end_pos,
        alt=alt,
    )


def _tra_events():
    return [
        _tra_event("seed", "chr1", 100, "chr2", 1000, "N[chr2:1000["),
        # Exact two-breakpoint tolerance boundary when delta=50 (TRA tolerance=100).
        _tra_event("boundary", "chr1", 200, "chr2", 1100, "N[chr2:1100["),
        _tra_event("outside", "chr1", 201, "chr2", 1101, "N[chr2:1101["),
        # Same physical adjacency represented from the reciprocal chromosome.
        _tra_event("reciprocal", "chr2", 1005, "chr1", 105, "]chr1:105]N"),
        # Coordinates are near but the known physical orientation is incompatible.
        _tra_event("incompatible", "chr1", 105, "chr2", 1005, "N]chr2:1005]"),
        # Unknown orientation must stay on the production legacy-coordinate fallback.
        _tra_event("symbolic", "chr1", 120, "chr2", 1020, "<TRA>"),
    ]


def _run_tra(*, force_full_scan: bool):
    merger = TRAMerger(delta=50, min_overlap_ratio=0.5, strand_consistency=True)
    if force_full_scan:
        # Returning no signature is the production-defined safe fallback that
        # causes merge_events() to scan every existing group while still using
        # the unchanged should_merge_tra() scientific rule.
        merger._adjacency_candidate_info = lambda event: (None, None)

    for event in _tra_events():
        merger.add_event(event)

    groups = []
    for pair_groups in merger.merge_events().values():
        groups.extend(pair_groups)
    return groups


def test_tra_candidate_pruning_preserves_full_scan_group_membership():
    optimized = _group_signature(_run_tra(force_full_scan=False))
    full_scan = _group_signature(_run_tra(force_full_scan=True))

    assert optimized == full_scan


def _bnd_event(event_id, chrom, pos, alt):
    return SimpleNamespace(
        sv_id=event_id,
        chrom=chrom,
        pos=pos,
        alt=alt,
    )


def _bnd_events():
    return [
        _bnd_event("seed", "chr1", 100, "N[chr2:1000["),
        # Exact local+remote delta boundary.
        _bnd_event("boundary", "chr1", 150, "N[chr2:1050["),
        _bnd_event("outside", "chr1", 151, "N[chr2:1051["),
        # Same coordinates, different BND pattern: must remain separate.
        _bnd_event("different.pattern", "chr1", 105, "N]chr2:1005]"),
        # Same numeric geometry, different target chromosome: must remain separate.
        _bnd_event("different.target", "chr1", 100, "N[chr3:1000["),
    ]


def _run_bnd(*, force_full_scan: bool):
    merger = BNDMerger(delta=50)
    if force_full_scan:
        # None is the production-defined full-scan fallback.  Scientific group
        # assignment still goes through should_merge_bnd().
        merger._bnd_candidate_key = lambda event: (None, None)

    for event in _bnd_events():
        merger.add_event(event)

    groups = []
    for pair_groups in merger.merge_events().values():
        groups.extend(pair_groups)
    return groups


def test_bnd_candidate_pruning_preserves_full_scan_group_membership():
    optimized = _group_signature(_run_bnd(force_full_scan=False))
    full_scan = _group_signature(_run_bnd(force_full_scan=True))

    assert optimized == full_scan
