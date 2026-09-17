"""Bridge caller-mode SVCF evidence into the shared sample-consensus model.

Caller-mode SVCF preserves every evidence block.  Consumers that need one
sample-level genotype (for example ``svcf2vcf`` or ``stat``) must therefore
reduce those evidence blocks without giving a caller multiple votes and
without reconstructing caller identity heuristically.

This module contains only the SVCF-to-consensus adapter.  The scientific state
machine remains in :mod:`octopusv.utils.sample_consensus`.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence

from octopusv.utils.sample_consensus import (
    SampleConsensusResult,
    resolve_sample_consensus,
)
from octopusv.utils.svcf_sample_parser import parse_svcf_sample_block


def _source_entries(info: Mapping[str, object]) -> list[str]:
    raw = info.get("SOURCES")
    if raw in (None, "", ".", True):
        return []
    return [entry.strip() for entry in str(raw).split(",")]


def resolve_caller_svcf_consensus(
    format_field: str,
    caller_segments: Sequence[str],
    info: Mapping[str, object],
) -> SampleConsensusResult:
    """Resolve one synthesized genotype from caller-mode SVCF evidence.

    ``SOURCES`` and evidence blocks are a positional contract.  When more than
    one evidence block is present, an exact caller identity is required for
    every block; a mismatch is an invalid/ambiguous caller-mode record and is
    rejected rather than guessed.

    A one-block record may omit ``SOURCES`` because there is no cross-caller
    ambiguity.  A private synthetic identity is used only to feed that single
    observation into the shared state machine.
    """

    segments = list(caller_segments)
    if not segments:
        return resolve_sample_consensus([])

    sources = _source_entries(info)

    if len(segments) == 1:
        if len(sources) > 1:
            raise ValueError(
                "Cannot synthesize caller genotype: SOURCES contains "
                f"{len(sources)} entries for 1 evidence block."
            )
        source = sources[0] if sources else "__single_evidence__"
        parsed = parse_svcf_sample_block(format_field, segments[0])
        return resolve_sample_consensus([(source, parsed.get("GT", "."))])

    if len(sources) != len(segments):
        raise ValueError(
            "Cannot synthesize caller genotype: SOURCES/evidence count "
            f"mismatch ({len(sources)} SOURCES entries vs "
            f"{len(segments)} evidence blocks)."
        )

    caller_genotypes = []
    for source, segment in zip(sources, segments):
        if source in ("", "."):
            raise ValueError(
                "Cannot synthesize caller genotype: one evidence block has "
                "no explicit caller identity in SOURCES."
            )
        parsed = parse_svcf_sample_block(format_field, segment)
        caller_genotypes.append((source, parsed.get("GT", ".")))

    return resolve_sample_consensus(caller_genotypes)
