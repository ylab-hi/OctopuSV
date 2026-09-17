"""Stable record ordering helpers for SVCF writers.

Output sorting is intentionally a representation-layer operation. It must not
change event grouping, representative selection, source/evidence binding, or
sample ordering.
"""

from __future__ import annotations

import re
from collections.abc import Iterable


_CONTIG_ID_RE = re.compile(r"^##contig=<ID=([^,>]+)")


def contig_order_from_meta_lines(meta_lines: Iterable[str]) -> list[str]:
    """Return contig IDs in the exact order declared by ``##contig`` lines."""
    contigs: list[str] = []
    seen: set[str] = set()

    for raw_line in meta_lines:
        match = _CONTIG_ID_RE.match(str(raw_line).strip())
        if match is None:
            continue
        contig = match.group(1)
        if contig not in seen:
            seen.add(contig)
            contigs.append(contig)

    return contigs


def sort_events_for_output(events, contig_order: Iterable[str]):
    """Return events in stable genomic output order.

    Declared header contig order is authoritative. Records on undeclared
    contigs follow declared contigs and are ordered by the exact contig string,
    without attempting aliases or natural-chromosome heuristics. Within a
    contig, records are ordered by POS, then END, SVTYPE, and record ID.
    The extra content-based tie-breakers make same-position output independent
    of upstream input order without changing any record contents.
    """
    rank: dict[str, int] = {}
    for index, contig in enumerate(contig_order):
        text = str(contig)
        if text not in rank:
            rank[text] = index

    def key(event):
        chrom = str(
            getattr(event, "chrom", None)
            or getattr(event, "start_chrom", "")
        )
        pos = getattr(event, "pos", None)
        if pos is None:
            pos = getattr(event, "start_pos", 0)
        pos = int(pos)

        info = getattr(event, "info", {}) or {}
        end_value = info.get("END", getattr(event, "end_pos", pos))
        try:
            end_key = (0, int(end_value))
        except (TypeError, ValueError):
            end_key = (1, str(end_value))

        svtype = str(info.get("SVTYPE", getattr(event, "sv_type", "")))
        record_id = str(
            getattr(event, "id", None)
            or getattr(event, "sv_id", "")
        )

        tie = (pos, end_key, svtype, record_id)
        if chrom in rank:
            return (0, rank[chrom], "", *tie)
        return (1, 0, chrom, *tie)

    return sorted(events, key=key)
