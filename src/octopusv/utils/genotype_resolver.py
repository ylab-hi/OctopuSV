"""Shared multi-caller genotype resolution.

In OctopuSV caller mode a single SVCF record may carry one FORMAT block per
supporting caller. When a single representative genotype is needed (genotype
statistics, or collapsing to a standard single-sample VCF column), the same
three-tier rule is applied so every consumer agrees:

  1. Majority vote across unique callers.
  2. Tie-break by AD variant-supporting reads (valid AD > missing AD).
  3. Tie-break by caller/source order (input file order).

If a merged SVCF record contains more than one evidence block from the same
source, only the first block from that source participates in genotype voting.
The additional blocks remain in the SVCF for traceability, but they do not give
one caller multiple votes.

This module is the single source of truth for that rule. Both
GenotypeAnalyzer (stat) and SVCFtoVCFConverter (svcf2vcf) import it.
"""

from collections import Counter


def extract_sources_order(info_field):
    """Extract SOURCES order (list of source names) from an INFO string."""
    try:
        for item in info_field.split(";"):
            if item.startswith("SOURCES="):
                return [s.strip() for s in item.split("=", 1)[1].split(",")]
        return []
    except Exception:
        return []


def unique_source_segments(caller_segments, info_field):
    """Return one evidence block per source, preserving source order.

    Current OctopuSV caller-mode SVCF records keep SOURCES and evidence blocks
    co-ordered. If one source contributes multiple records to a merged event,
    the additional evidence blocks are preserved in the SVCF, but a caller must
    not receive multiple genotype votes merely because it contributed multiple
    nearby records.

    For a well-formed record, keep the first evidence block for each distinct
    SOURCES entry. If SOURCES is absent or its length does not match the number
    of evidence blocks, return all blocks unchanged rather than guessing.

    Returns:
        List of ``(original_index, segment)`` pairs.
    """
    sources = extract_sources_order(info_field)

    if not sources or len(sources) != len(caller_segments):
        return list(enumerate(caller_segments))

    seen = set()
    selected = []

    for index, (source, segment) in enumerate(zip(sources, caller_segments)):
        if source in seen:
            continue

        seen.add(source)
        selected.append((index, segment))

    return selected


def extract_variant_support(ad_value):
    """Variant-supporting read count from an AD field (its second number).

    Returns -1 for missing/empty values so they are deprioritized.
    """
    try:
        if not ad_value or ad_value in (".", ""):
            return -1
        if "," in ad_value:
            parts = ad_value.split(",")
            if len(parts) >= 2:
                second = parts[1].strip()
                if second in (".", ""):
                    return -1
                return int(second)
            first = parts[0].strip()
            if first in (".", ""):
                return -1
            return int(first)
        if ad_value not in (".", ""):
            return int(ad_value)
        return -1
    except (ValueError, IndexError):
        return -1


def is_better_ad_support(new_ad, current_ad):
    """Valid AD (>=0) beats missing AD (-1); among valid, higher wins."""
    if current_ad < 0 and new_ad >= 0:
        return True
    if new_ad < 0 and current_ad >= 0:
        return False
    return new_ad > current_ad


def resolve_multi_caller_genotype(format_field, caller_segments, info_field):
    """Resolve a single representative genotype from caller evidence blocks.

    Each unique source contributes at most one vote. If the same source appears
    multiple times, the first block for that source is used for genotype
    resolution so preserving extra evidence does not change caller weighting.

    Args:
        format_field: the FORMAT string, e.g. "GT:AD:LN:ST:QV:TY:ID:SC:REF:ALT:CO".
        caller_segments: raw caller evidence strings.
        info_field: the record INFO string containing SOURCES when available.

    Returns:
        The winning genotype string (e.g. "1/1"), or None if GT is absent.
    """
    format_keys = format_field.split(":")
    if "GT" not in format_keys:
        return None

    gt_index = format_keys.index("GT")
    ad_index = format_keys.index("AD") if "AD" in format_keys else None

    selected_segments = unique_source_segments(caller_segments, info_field)

    genotype_data = []  # (gt, ad_support, original_caller_index)
    for original_index, seg in selected_segments:
        fields = seg.split(":")
        if gt_index < len(fields):
            gt = fields[gt_index]
            ad_support = 0
            if ad_index is not None and ad_index < len(fields):
                ad_support = extract_variant_support(fields[ad_index])
            genotype_data.append((gt, ad_support, original_index))

    if not genotype_data:
        return None

    # Tier 1: majority vote across unique sources.
    vote_counts = Counter(item[0] for item in genotype_data)
    max_votes = max(vote_counts.values())
    tied = [gt for gt, count in vote_counts.items() if count == max_votes]
    if len(tied) == 1:
        return tied[0]

    # Tier 2: AD support among tied genotypes.
    tied_ad = {}
    for gt, ad_support, index in genotype_data:
        if gt in tied:
            if gt not in tied_ad or is_better_ad_support(
                ad_support,
                tied_ad[gt][0],
            ):
                tied_ad[gt] = (ad_support, index)

    valid_ad = {gt: value for gt, value in tied_ad.items() if value[0] >= 0}
    if valid_ad:
        max_ad = max(ad for ad, _ in valid_ad.values())
        ad_winners = [
            gt
            for gt, (ad, _) in valid_ad.items()
            if ad == max_ad
        ]
    else:
        ad_winners = list(tied_ad.keys())

    if len(ad_winners) == 1:
        return ad_winners[0]

    # Tier 3: earliest source/evidence block among remaining tied genotypes.
    earliest = min(tied_ad[gt][1] for gt in ad_winners)
    for gt in ad_winners:
        if tied_ad[gt][1] == earliest:
            return gt

    return ad_winners[0]
