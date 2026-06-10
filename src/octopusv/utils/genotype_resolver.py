"""Shared multi-caller genotype resolution.

In OctopuSV caller mode a single SVCF record may carry one FORMAT block per
supporting caller. When a single representative genotype is needed (genotype
statistics, or collapsing to a standard single-sample VCF column), the same
three-tier rule is applied so every consumer agrees:

  1. Majority vote across caller genotypes.
  2. Tie-break by AD variant-supporting reads (valid AD > missing AD).
  3. Tie-break by SOURCES order (input file order); then first caller block.

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
    """Resolve a single representative genotype from multiple caller blocks.

    Args:
        format_field: the FORMAT string, e.g. "GT:AD:LN:ST:QV:TY:ID:SC:REF:ALT:CO".
        caller_segments: list of raw sample-column strings, one per caller.
        info_field: the record INFO string (used for SOURCES-order tie-break).

    Returns:
        The winning genotype string (e.g. "1/1"), or None if GT is absent.
    """
    format_keys = format_field.split(":")
    if "GT" not in format_keys:
        return None

    gt_index = format_keys.index("GT")
    ad_index = format_keys.index("AD") if "AD" in format_keys else None

    genotype_data = []  # (gt, ad_support, caller_index)
    for i, seg in enumerate(caller_segments):
        fields = seg.split(":")
        if gt_index < len(fields):
            gt = fields[gt_index]
            ad_support = 0
            if ad_index is not None and ad_index < len(fields):
                ad_support = extract_variant_support(fields[ad_index])
            genotype_data.append((gt, ad_support, i))

    if not genotype_data:
        return None

    # Tier 1: majority vote.
    vote_counts = Counter(item[0] for item in genotype_data)
    max_votes = max(vote_counts.values())
    tied = [gt for gt, c in vote_counts.items() if c == max_votes]
    if len(tied) == 1:
        return tied[0]

    # Tier 2: AD support among tied genotypes.
    tied_ad = {}
    for gt, ad_support, idx in genotype_data:
        if gt in tied:
            if gt not in tied_ad or is_better_ad_support(ad_support, tied_ad[gt][0]):
                tied_ad[gt] = (ad_support, idx)

    valid_ad = {gt: v for gt, v in tied_ad.items() if v[0] >= 0}
    if valid_ad:
        max_ad = max(ad for ad, _ in valid_ad.values())
        ad_winners = [gt for gt, (ad, _) in valid_ad.items() if ad == max_ad]
    else:
        ad_winners = list(tied_ad.keys())

    if len(ad_winners) == 1:
        return ad_winners[0]

    # Tier 3: SOURCES order (input file order) -> earliest caller block.
    sources_order = extract_sources_order(info_field)
    if sources_order:
        return min(ad_winners, key=lambda gt: tied_ad[gt][1])

    # Fallback: earliest caller block among remaining.
    earliest = min(tied_ad[gt][1] for gt in ad_winners)
    for gt in ad_winners:
        if tied_ad[gt][1] == earliest:
            return gt
    return ad_winners[0]
