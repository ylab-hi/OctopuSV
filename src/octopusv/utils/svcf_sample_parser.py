"""Shared parser for SVCF FORMAT evidence blocks.

SVCF uses the fixed evidence tail::

    ID:SC:REF:ALT:CO

Some values in that tail may legally contain ``:``:

* source record IDs (for example Manta IDs)
* BND/TRA ALT strings (for example ``N]chr2:12345]``)
* symbolic ALT subtypes (for example ``<INS:ME:ALU>``)

A plain ``str.split(\":\")`` therefore cannot recover the fields reliably.
This module provides the single parser used by SVCF readers that need the
full FORMAT/evidence structure.
"""

from __future__ import annotations

import re
from collections.abc import Sequence


STANDARD_EVIDENCE_TAIL = ("ID", "SC", "REF", "ALT", "CO")

# VCF breakend ALT forms have one remote ``chrom:pos`` coordinate enclosed
# by matching '[' or ']' brackets.  The local replacement sequence may use
# more than A/C/G/T/N, so do not unnecessarily restrict its alphabet here.
_BND_ALT_PATTERN = re.compile(
    r"^[^:\[\]]*([\[\]])[^:\[\]]+:\d+\1[^:\[\]]*$"
)

# Symbolic ALT identifiers may contain colon-delimited subtypes, e.g.
# <INS:ME:ALU>.  SVCF preserves the original source ALT value.
_SYMBOLIC_ALT_PATTERN = re.compile(r"^<[^<>]+>$")


def _format_keys(format_field_or_keys: str | Sequence[str]) -> list[str]:
    if isinstance(format_field_or_keys, str):
        if not format_field_or_keys:
            return []
        return format_field_or_keys.split(":")

    return [str(key) for key in format_field_or_keys]


def _generic_parse(keys: list[str], parts: list[str]) -> dict[str, str]:
    """Conservative parser for non-standard FORMAT layouts.

    The final FORMAT key absorbs overflow so a colon in the final value does
    not silently create extra fields. Missing values are filled with ``.``.
    """
    result: dict[str, str] = {}

    for index, key in enumerate(keys):
        if index >= len(parts):
            result[key] = "."
            continue

        if index == len(keys) - 1:
            result[key] = ":".join(parts[index:])
        else:
            result[key] = parts[index]

    return result


def _looks_like_bnd_alt(value: str) -> bool:
    return bool(_BND_ALT_PATTERN.fullmatch(value))


def _symbolic_alt_suffix_start(parts: list[str]) -> int | None:
    """Return the start index of a symbolic ALT occupying the list suffix."""
    if not parts or not parts[-1].endswith(">"):
        return None

    # Walk from the right because ALT is the final field before CO.  The first
    # suffix that forms one complete symbolic allele is the ALT value.
    for index in range(len(parts) - 1, -1, -1):
        if not parts[index].startswith("<"):
            continue

        candidate = ":".join(parts[index:])
        if _SYMBOLIC_ALT_PATTERN.fullmatch(candidate):
            return index

    return None


def parse_svcf_sample_block(
    format_field_or_keys: str | Sequence[str],
    block: str | None,
) -> dict[str, str]:
    """Parse one SVCF sample/evidence block.

    For OctopuSV's fixed SVCF schema, fields before ``ID`` are parsed from the
    left and ``ID:SC:REF:ALT:CO`` is resolved from both ends. This preserves
    colon-containing source IDs, BND/TRA ALT strings, and colon-containing
    symbolic ALT subtypes without changing the SVCF format.

    Args:
        format_field_or_keys: FORMAT string or an ordered list of FORMAT keys.
        block: One tab-delimited SVCF evidence/sample value.

    Returns:
        Dictionary containing one value for every FORMAT key. Missing values
        are represented by ``.``. ``{}`` is returned when FORMAT is empty or
        the block is ``None``.
    """
    keys = _format_keys(format_field_or_keys)

    if not keys or block is None:
        return {}

    parts = str(block).split(":")

    if tuple(keys[-5:]) != STANDARD_EVIDENCE_TAIL:
        return _generic_parse(keys, parts)

    head_count = len(keys) - len(STANDARD_EVIDENCE_TAIL)

    # A short block cannot contain an expanded colon-bearing ID/ALT. Fill it
    # positionally rather than inventing field boundaries.
    if len(parts) < len(keys):
        return {
            key: parts[index] if index < len(parts) else "."
            for index, key in enumerate(keys)
        }

    result: dict[str, str] = {}

    for index in range(head_count):
        result[keys[index]] = parts[index]

    tail_parts = parts[head_count:]

    if len(tail_parts) < 5:
        # Defensive fallback; the length check above should normally prevent
        # this path for the standard SVCF schema.
        for index, key in enumerate(STANDARD_EVIDENCE_TAIL):
            result[key] = tail_parts[index] if index < len(tail_parts) else "."
        return result

    co_value = tail_parts[-1]
    before_co = tail_parts[:-1]

    # ALT is the suffix immediately before CO.  A BND ALT occupies exactly two
    # colon tokens because of its remote chrom:pos coordinate.  A symbolic ALT
    # such as <INS:ME:ALU> may occupy two or more tokens.  Ordinary sequence or
    # simple symbolic ALT values occupy one token.
    alt_start = len(before_co) - 1
    alt_value = before_co[-1]

    if len(before_co) >= 2:
        bnd_candidate = f"{before_co[-2]}:{before_co[-1]}"
        if _looks_like_bnd_alt(bnd_candidate):
            alt_start = len(before_co) - 2
            alt_value = bnd_candidate
        else:
            symbolic_start = _symbolic_alt_suffix_start(before_co)
            if symbolic_start is not None:
                alt_start = symbolic_start
                alt_value = ":".join(before_co[symbolic_start:])

    ref_index = alt_start - 1
    sc_index = ref_index - 1

    if sc_index < 1:
        # There is not enough structure left for ID + SC + REF + ALT. Do not
        # guess; preserve a deterministic positional interpretation.
        return _generic_parse(keys, parts)

    id_tokens = before_co[:sc_index]

    result["ID"] = ":".join(id_tokens) if id_tokens else "."
    result["SC"] = before_co[sc_index]
    result["REF"] = before_co[ref_index]
    result["ALT"] = alt_value
    result["CO"] = co_value

    return result
