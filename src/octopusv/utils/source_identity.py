"""Shared helpers for explicit SVCF source identities.

Source labels are factual identities.  OctopuSV may report likely case-only
CLI typos, but it must never normalize, infer, or rewrite those identities.
"""

from __future__ import annotations

from octopusv.utils.svcf_sample_parser import parse_svcf_sample_block


MISSING_SOURCE_VALUES = {None, "", ".", True}


def explicit_record_sources(info: dict, fields: list[str]) -> tuple[set[str], str]:
    """Return explicit record-level sources and where they came from.

    Merged records use INFO/SOURCES.  A single-evidence caller record may
    instead carry source identity in FORMAT/SC.  No filename/ID/header-label
    guessing is performed.
    """
    sources_value = info.get("SOURCES")
    if sources_value not in MISSING_SOURCE_VALUES:
        sources = {
            item.strip()
            for item in str(sources_value).split(",")
            if item.strip() and item.strip() != "."
        }
        if sources:
            return sources, "INFO/SOURCES"

    if len(fields) == 10:
        parsed = parse_svcf_sample_block(fields[8], fields[9])
        source = parsed.get("SC")
        if source not in (None, "", ".", "unknown"):
            return {str(source)}, "FORMAT/SC"

    raise ValueError(
        "Source-based operations require explicit source identity in "
        "INFO/SOURCES or, for a single-evidence record, FORMAT/SC. "
        "OctopuSV does not infer source identity from record IDs, file "
        "names, or #CHROM labels."
    )


def case_only_source_candidates(requested: str, observed_sources) -> list[str]:
    """Return observed labels differing from ``requested`` only by case."""
    return sorted(
        observed
        for observed in set(observed_sources)
        if observed != requested and observed.casefold() == requested.casefold()
    )


def format_case_only_source_error(requested: str, candidates: list[str]) -> str:
    """Return the shared actionable error for a case-only source typo."""
    if len(candidates) == 1:
        suggestion = f"Did you mean {candidates[0]!r}?"
    else:
        suggestion = "Case-sensitive candidates are: " + ", ".join(
            repr(candidate) for candidate in candidates
        )
    return (
        f"No exact source {requested!r}. {suggestion} "
        "Source identity is case-sensitive."
    )
