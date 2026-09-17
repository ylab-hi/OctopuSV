"""Coordinate parsing helpers for SVCF evidence fields.

The SVCF ``CO`` value has the form::

    startChrom_startPos-endChrom_endPos

Contig names may themselves contain underscores or hyphens.  Parsing must
therefore not assume that the first underscore or first hyphen is a delimiter.
This module keeps the parser conservative: it returns a coordinate tuple only
when the text has exactly one structurally valid interpretation.
"""

from __future__ import annotations


def _parse_endpoint(value: str) -> tuple[str, int] | None:
    """Parse ``chrom_pos`` by splitting at the rightmost underscore."""
    chrom, sep, pos_text = value.rpartition("_")
    if not sep or not chrom or not pos_text.isdigit():
        return None

    pos = int(pos_text)
    if pos < 0:
        return None

    return chrom, pos


def parse_svcf_co(value: str | None) -> tuple[str, int, str, int] | None:
    """Parse one SVCF ``CO`` value.

    Coordinate positions are parsed structurally and may be zero here; fixed-field
    POS validity is checked separately by the SVCF validator.

    Hyphens are allowed inside contig names, so every hyphen is considered as
    a possible start/end separator.  A result is returned only when exactly one
    split yields two valid ``chrom_pos`` endpoints.  Ambiguous text is rejected
    rather than guessed.
    """
    if value in (None, "", "."):
        return None

    text = str(value)
    candidates: list[tuple[str, int, str, int]] = []

    for index, char in enumerate(text):
        if char != "-":
            continue

        left = _parse_endpoint(text[:index])
        right = _parse_endpoint(text[index + 1 :])
        if left is None or right is None:
            continue

        candidate = (left[0], left[1], right[0], right[1])
        if candidate not in candidates:
            candidates.append(candidate)

    if len(candidates) != 1:
        return None

    return candidates[0]
