"""Executable contract for OctopuSV 1.0 sample-level genotype consensus.

These tests intentionally define the sample-consensus semantics before the
production resolver is implemented.  They are pure state-machine tests: no
files, merge engine, writer, or converter is involved.

Expected production API (to be added later):

    octopusv.utils.sample_consensus.reduce_caller_genotypes(genotypes)
    octopusv.utils.sample_consensus.resolve_sample_consensus(caller_genotypes)

``caller_genotypes`` is an iterable of ``(source, GT)`` pairs.  Multiple pairs
with the same source are multiple evidence records from one caller and must be
reduced to one caller state before cross-caller consensus.
"""

from __future__ import annotations

from itertools import permutations

import pytest


try:
    from octopusv.utils.sample_consensus import (
        reduce_caller_genotypes,
        resolve_sample_consensus,
    )
except ModuleNotFoundError:
    reduce_caller_genotypes = None
    resolve_sample_consensus = None


NO_VOTE = "NO_VOTE"
ABSENT = "ABSENT"
CARRIER_HET = "CARRIER_HET"
CARRIER_HOM = "CARRIER_HOM"
CARRIER_UNKNOWN = "CARRIER_UNKNOWN"


def _require_api():
    assert reduce_caller_genotypes is not None, (
        "OctopuSV 1.0 sample-consensus API is not implemented yet: "
        "octopusv.utils.sample_consensus.reduce_caller_genotypes"
    )
    assert resolve_sample_consensus is not None, (
        "OctopuSV 1.0 sample-consensus API is not implemented yet: "
        "octopusv.utils.sample_consensus.resolve_sample_consensus"
    )


def _state(genotypes):
    _require_api()
    state = reduce_caller_genotypes(genotypes)
    # Permit either a plain string or a string-like Enum while keeping the
    # scientific contract independent of implementation style.
    return getattr(state, "value", state)


def _resolve(caller_genotypes):
    _require_api()
    result = resolve_sample_consensus(caller_genotypes)

    # Keep the public contract intentionally small: the resolver must expose
    # the synthesized GT plus the two cohort-QC counts planned for SVCF 1.1.
    assert hasattr(result, "gt")
    assert hasattr(result, "unique_carrier_callers")
    assert hasattr(result, "valid_caller_votes")
    return result


@pytest.mark.parametrize(
    ("genotypes", "expected"),
    [
        (["./."], NO_VOTE),
        (["."], NO_VOTE),
        (["0/0"], ABSENT),
        (["0|0"], ABSENT),
        (["0/1"], CARRIER_HET),
        (["1/0"], CARRIER_HET),
        (["0|1"], CARRIER_HET),
        (["1|0"], CARRIER_HET),
        (["1/1"], CARRIER_HOM),
        (["1|1"], CARRIER_HOM),
    ],
)
def test_single_evidence_reduces_to_expected_caller_state(genotypes, expected):
    assert _state(genotypes) == expected


def test_missing_evidence_does_not_override_a_valid_caller_call():
    assert _state(["./.", "0/1", "."]) == CARRIER_HET


def test_same_caller_repeated_het_evidence_is_one_het_state():
    assert _state(["0/1", "1/0", "0|1"]) == CARRIER_HET


def test_same_caller_repeated_hom_alt_evidence_is_one_hom_state():
    assert _state(["1/1", "1|1"]) == CARRIER_HOM


def test_same_caller_het_hom_disagreement_preserves_carrier_presence():
    assert _state(["0/1", "1/1"]) == CARRIER_UNKNOWN


def test_same_caller_absent_and_carrier_evidence_is_no_vote():
    assert _state(["0/0", "0/1"]) == NO_VOTE


def test_no_valid_caller_votes_returns_missing():
    result = _resolve([
        ("callerA", "./."),
        ("callerB", "."),
    ])

    assert result.gt == "./."
    assert result.unique_carrier_callers == 0
    assert result.valid_caller_votes == 0


def test_one_carrier_vote_is_not_diluted_by_missing_callers():
    result = _resolve([
        ("callerA", "0/1"),
        ("callerB", "./."),
        ("callerC", "."),
    ])

    assert result.gt == "0/1"
    assert result.unique_carrier_callers == 1
    assert result.valid_caller_votes == 1


def test_all_absent_votes_return_hom_reference():
    result = _resolve([
        ("callerA", "0/0"),
        ("callerB", "0|0"),
    ])

    assert result.gt == "0/0"
    assert result.unique_carrier_callers == 0
    assert result.valid_caller_votes == 2


def test_all_het_carriers_return_het():
    result = _resolve([
        ("callerA", "0/1"),
        ("callerB", "1/0"),
        ("callerC", "0|1"),
    ])

    assert result.gt == "0/1"
    assert result.unique_carrier_callers == 3
    assert result.valid_caller_votes == 3


def test_all_hom_alt_carriers_return_hom_alt():
    result = _resolve([
        ("callerA", "1/1"),
        ("callerB", "1|1"),
    ])

    assert result.gt == "1/1"
    assert result.unique_carrier_callers == 2
    assert result.valid_caller_votes == 2


def test_het_hom_disagreement_returns_half_call_not_false_het_or_missing():
    result = _resolve([
        ("callerA", "0/1"),
        ("callerB", "1/1"),
    ])

    assert result.gt == "1/."
    assert result.unique_carrier_callers == 2
    assert result.valid_caller_votes == 2


def test_unknown_carrier_zygosity_returns_half_call():
    result = _resolve([
        ("callerA", "0/1"),
        ("callerA", "1/1"),
        ("callerB", "0/1"),
    ])

    assert result.gt == "1/."
    assert result.unique_carrier_callers == 2
    assert result.valid_caller_votes == 2


def test_carrier_presence_strict_majority_wins():
    result = _resolve([
        ("callerA", "0/1"),
        ("callerB", "1/1"),
        ("callerC", "0/0"),
    ])

    # Two unique callers support presence, one supports absence.  Presence wins;
    # HET/HOM disagreement then makes zygosity partially unknown.
    assert result.gt == "1/."
    assert result.unique_carrier_callers == 2
    assert result.valid_caller_votes == 3


def test_absent_presence_strict_majority_wins():
    result = _resolve([
        ("callerA", "0/0"),
        ("callerB", "0/0"),
        ("callerC", "1/1"),
    ])

    assert result.gt == "0/0"
    assert result.unique_carrier_callers == 1
    assert result.valid_caller_votes == 3


def test_true_presence_tie_returns_missing():
    result = _resolve([
        ("callerA", "0/1"),
        ("callerB", "0/0"),
    ])

    assert result.gt == "./."
    assert result.unique_carrier_callers == 1
    assert result.valid_caller_votes == 2


def test_same_caller_multiple_evidence_contributes_exactly_one_caller_state():
    result = _resolve([
        ("callerA", "0/1"),
        ("callerA", "1/1"),
        ("callerA", "0/1"),
        ("callerB", "0/1"),
    ])

    # callerA is one CARRIER_UNKNOWN state, not three votes.
    assert result.gt == "1/."
    assert result.unique_carrier_callers == 2
    assert result.valid_caller_votes == 2


def test_internally_contradictory_caller_is_no_vote_not_absent():
    result = _resolve([
        ("callerA", "0/0"),
        ("callerA", "0/1"),
        ("callerB", "0/1"),
    ])

    # callerA becomes NO_VOTE. callerB alone establishes carrier presence.
    assert result.gt == "0/1"
    assert result.unique_carrier_callers == 1
    assert result.valid_caller_votes == 1


def test_consensus_is_independent_of_caller_and_evidence_order():
    evidence = [
        ("callerA", "0/1"),
        ("callerA", "1/1"),
        ("callerB", "0/1"),
        ("callerC", "0/0"),
    ]

    expected = None
    for permuted in permutations(evidence):
        result = _resolve(permuted)
        observed = (
            result.gt,
            result.unique_carrier_callers,
            result.valid_caller_votes,
        )
        if expected is None:
            expected = observed
        else:
            assert observed == expected

    assert expected == ("1/.", 2, 3)


def test_haploid_reference_and_alt_are_supported_without_silent_no_vote():
    assert _state(["0"]) == ABSENT
    assert _state(["1"]) == CARRIER_UNKNOWN

    ref_result = _resolve([
        ("callerA", "0"),
        ("callerB", "0"),
    ])
    assert ref_result.gt == "0"
    assert ref_result.unique_carrier_callers == 0
    assert ref_result.valid_caller_votes == 2

    alt_result = _resolve([
        ("callerA", "1"),
        ("callerB", "1"),
    ])
    assert alt_result.gt == "1"
    assert alt_result.unique_carrier_callers == 2
    assert alt_result.valid_caller_votes == 2


def test_mixed_haploid_and_diploid_carrier_evidence_preserves_presence_only():
    result = _resolve([
        ("callerA", "1"),
        ("callerB", "0/1"),
    ])

    assert result.gt == "1/."
    assert result.unique_carrier_callers == 2
    assert result.valid_caller_votes == 2


def test_multiallelic_or_unsupported_genotypes_do_not_vote():
    assert _state(["0/2"]) == NO_VOTE
    assert _state(["1/2"]) == NO_VOTE
    assert _state(["2"]) == NO_VOTE

    result = _resolve([
        ("callerA", "0/2"),
        ("callerB", "1/2"),
    ])
    assert result.gt == "./."
    assert result.unique_carrier_callers == 0
    assert result.valid_caller_votes == 0


def test_explicit_zero_zero_is_a_real_absence_vote_not_a_layout_placeholder():
    # This test intentionally locks the resolver-side contract: if an upstream
    # layer passes 0/0 here, it means a caller explicitly asserted absence.
    # Sample-mode fixed-width placeholder columns therefore MUST be filtered
    # out before this resolver is called.
    result = _resolve([
        ("callerA", "0/0"),
        ("callerB", "0/0"),
        ("callerC", "0/1"),
    ])

    assert result.gt == "0/0"
    assert result.unique_carrier_callers == 1
    assert result.valid_caller_votes == 3


def test_haploid_consensus_is_independent_of_caller_order():
    evidence = [
        ("callerA", "1"),
        ("callerB", "1"),
        ("callerC", "0"),
    ]

    observed = {
        (
            _resolve(permuted).gt,
            _resolve(permuted).unique_carrier_callers,
            _resolve(permuted).valid_caller_votes,
        )
        for permuted in permutations(evidence)
    }
    assert observed == {("1", 2, 3)}
