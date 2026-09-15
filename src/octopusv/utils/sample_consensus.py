"""Order-independent sample-level genotype consensus for OctopuSV 1.0.

Caller-mode SVCF preserves every caller evidence block. Sample mode instead
needs one synthesized genotype per biological sample. This module implements
that synthesis without using caller input order, caller-specific AD values, or
other heuristic tie-breakers.

The resolver deliberately separates two questions:

1. Does this caller support carrier presence, absence, or no reliable vote?
2. If carrier presence is established across callers, is zygosity resolved?

Each unique caller contributes exactly one state. Missing callers are not
interpreted as ``0/0``. If carrier presence wins but diploid carrier callers
disagree on heterozygous versus homozygous-alt status, the synthesized GT is
``1/.``: at least one ALT allele is established while the second allele is
unknown.

Input contract
--------------
``resolve_sample_consensus`` accepts *caller-mode evidence only*: every
``(source, GT)`` pair must correspond to a real evidence block emitted by a
caller. Sample-mode layout placeholders must never be passed to this resolver.
In particular, if ``0/0`` is supplied it is interpreted as an explicit caller
claim of reference/absence, not as a fixed-width column placeholder. A caller
that did not contribute evidence is simply absent from the input (or has only a
missing GT and therefore contributes ``NO_VOTE``).

The genotype vocabulary is biallelic. Haploid ``0`` and ``1`` are supported;
multiallelic genotypes such as ``0/2`` or ``1/2`` are intentionally ``NO_VOTE``
rather than being guessed into the biallelic model.
"""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum
from typing import Hashable, Iterable, List, Tuple


class CallerState(str, Enum):
    """One genotype state contributed by one unique caller."""

    NO_VOTE = "NO_VOTE"
    ABSENT = "ABSENT"
    CARRIER_HET = "CARRIER_HET"
    CARRIER_HOM = "CARRIER_HOM"
    CARRIER_UNKNOWN = "CARRIER_UNKNOWN"


@dataclass(frozen=True)
class SampleConsensusResult:
    """Synthesized sample-level genotype plus minimal cohort-QC counts."""

    gt: str
    unique_carrier_callers: int
    valid_caller_votes: int


_MISSING_GT = {"", ".", "./.", ".|."}
_HAP_REF_GT = {"0"}
_HAP_ALT_GT = {"1"}
_HET_GT = {"0/1", "1/0", "0|1", "1|0"}
_HOM_REF_GT = {"0/0", "0|0"}
_HOM_ALT_GT = {"1/1", "1|1"}
_HALF_CARRIER_GT = {"1/.", "./1", "1|.", ".|1"}


def _normalized_gt(genotype: object) -> str:
    if genotype is None:
        return ""
    return str(genotype).strip()


def _classify_single_genotype(genotype: object) -> CallerState:
    """Classify one supported biallelic GT into the caller-state vocabulary.

    Haploid ``0`` is explicit absence and haploid ``1`` establishes carrier
    presence with no diploid zygosity claim, so it maps to
    ``CARRIER_UNKNOWN`` at the state layer. Unsupported or malformed GTs,
    including multiallelic calls, contribute ``NO_VOTE``.
    """

    gt = _normalized_gt(genotype)
    if gt in _MISSING_GT:
        return CallerState.NO_VOTE
    if gt in _HAP_REF_GT or gt in _HOM_REF_GT:
        return CallerState.ABSENT
    if gt in _HAP_ALT_GT:
        return CallerState.CARRIER_UNKNOWN
    if gt in _HET_GT:
        return CallerState.CARRIER_HET
    if gt in _HOM_ALT_GT:
        return CallerState.CARRIER_HOM
    if gt in _HALF_CARRIER_GT:
        return CallerState.CARRIER_UNKNOWN

    return CallerState.NO_VOTE


def _informative_supported_gts(genotypes: Iterable[object]) -> List[str]:
    """Return supported, non-missing GT strings used by one caller."""

    supported = (
        _HAP_REF_GT
        | _HAP_ALT_GT
        | _HET_GT
        | _HOM_REF_GT
        | _HOM_ALT_GT
        | _HALF_CARRIER_GT
    )
    result: List[str] = []
    for genotype in genotypes:
        gt = _normalized_gt(genotype)
        if gt in supported:
            result.append(gt)
    return result


def _is_haploid_alt_only(genotypes: Iterable[object]) -> bool:
    informative = _informative_supported_gts(genotypes)
    return bool(informative) and all(gt in _HAP_ALT_GT for gt in informative)


def _is_haploid_ref_only(genotypes: Iterable[object]) -> bool:
    informative = _informative_supported_gts(genotypes)
    return bool(informative) and all(gt in _HAP_REF_GT for gt in informative)


def reduce_caller_genotypes(genotypes: Iterable[object]) -> CallerState:
    """Reduce all evidence GTs from one caller to exactly one caller state.

    Missing or unsupported evidence is ignored when at least one informative
    supported genotype exists. Multiple non-reference evidence blocks preserve
    carrier presence even when they disagree on zygosity. A caller that
    contains both explicit absence (haploid ``0`` or diploid ``0/0``) and
    non-reference evidence is internally contradictory and therefore
    contributes ``NO_VOTE``.
    """

    informative: List[CallerState] = []
    for genotype in genotypes:
        state = _classify_single_genotype(genotype)
        if state is not CallerState.NO_VOTE:
            informative.append(state)

    if not informative:
        return CallerState.NO_VOTE

    states = set(informative)

    carrier_states = {
        CallerState.CARRIER_HET,
        CallerState.CARRIER_HOM,
        CallerState.CARRIER_UNKNOWN,
    }
    if CallerState.ABSENT in states and states.intersection(carrier_states):
        return CallerState.NO_VOTE

    if states == {CallerState.ABSENT}:
        return CallerState.ABSENT
    if states == {CallerState.CARRIER_HET}:
        return CallerState.CARRIER_HET
    if states == {CallerState.CARRIER_HOM}:
        return CallerState.CARRIER_HOM

    if states.issubset(carrier_states):
        return CallerState.CARRIER_UNKNOWN

    return CallerState.NO_VOTE


def resolve_sample_consensus(
    caller_genotypes: Iterable[Tuple[Hashable, object]],
) -> SampleConsensusResult:
    """Resolve an order-independent sample-level genotype across callers.

    Args:
        caller_genotypes: Iterable of ``(source, GT)`` pairs from caller-mode
            evidence blocks. Multiple pairs with the same source are evidence
            records from the same caller and are reduced to one caller state
            before cross-caller synthesis. Fixed-width sample-mode placeholders
            are outside this API contract and must be excluded upstream.

    Returns:
        ``SampleConsensusResult`` containing the synthesized GT, the number of
        unique callers supporting carrier presence, and the number of unique
        callers contributing a valid presence vote.

    Rules:
        * ``NO_VOTE`` callers are excluded from presence voting.
        * Carrier presence versus explicit absence is decided by strict
          majority among valid caller states; a true tie is ``./.``.
        * When carrier presence wins, all diploid HET -> ``0/1`` and all
          diploid HOM-alt -> ``1/1``. Diploid zygosity disagreement or
          uncertainty -> ``1/.``.
        * When every winning carrier vote is explicitly haploid ALT, the
          synthesized GT remains haploid ``1``. Likewise, an all-haploid
          reference/absence result remains ``0``. Mixed haploid/diploid carrier
          evidence preserves carrier status but reports unresolved diploid
          zygosity as ``1/.``.
    """

    grouped: dict[Hashable, List[object]] = {}
    for source, genotype in caller_genotypes:
        grouped.setdefault(source, []).append(genotype)

    caller_entries = [
        (reduce_caller_genotypes(genotypes), genotypes)
        for genotypes in grouped.values()
    ]

    carrier_states = {
        CallerState.CARRIER_HET,
        CallerState.CARRIER_HOM,
        CallerState.CARRIER_UNKNOWN,
    }

    carrier_votes = [
        (state, genotypes)
        for state, genotypes in caller_entries
        if state in carrier_states
    ]
    absent_votes = [
        (state, genotypes)
        for state, genotypes in caller_entries
        if state is CallerState.ABSENT
    ]

    unique_carrier_callers = len(carrier_votes)
    valid_caller_votes = unique_carrier_callers + len(absent_votes)

    if valid_caller_votes == 0:
        return SampleConsensusResult("./.", 0, 0)

    if unique_carrier_callers == len(absent_votes):
        return SampleConsensusResult(
            "./.",
            unique_carrier_callers,
            valid_caller_votes,
        )

    if unique_carrier_callers < len(absent_votes):
        # Preserve haploid reference output only when every caller that actually
        # supports absence does so with an explicitly haploid reference GT.
        gt = (
            "0"
            if all(_is_haploid_ref_only(gts) for _, gts in absent_votes)
            else "0/0"
        )
        return SampleConsensusResult(
            gt,
            unique_carrier_callers,
            valid_caller_votes,
        )

    # Carrier presence has strict majority. Preserve haploid ALT output only
    # when every carrier voter is explicitly haploid. Otherwise resolve the
    # diploid zygosity claim conservatively.
    if all(_is_haploid_alt_only(gts) for _, gts in carrier_votes):
        gt = "1"
    else:
        carrier_state_set = {state for state, _ in carrier_votes}
        if carrier_state_set == {CallerState.CARRIER_HET}:
            gt = "0/1"
        elif carrier_state_set == {CallerState.CARRIER_HOM}:
            gt = "1/1"
        else:
            gt = "1/."

    return SampleConsensusResult(
        gt,
        unique_carrier_callers,
        valid_caller_votes,
    )
