"""Shared downstream semantics for synthesized sample-mode SVCF calls.

SVCF 1.1 sample mode uses fixed-width sample columns.  When a biological
sample does not contribute a merged event, ``MultiSampleWriter`` writes a
layout placeholder with ``GT=0/0``, ``UC=0``, ``UV=0`` and no source/event
identity (``SC=.`` / ``ID=.``).  That placeholder is deliberately distinct
from both:

* an evidence-backed homozygous-reference call (for example ``GT=0/0, UV=2``),
* an unresolved synthesized call that *does* have evidence but no valid vote
  (for example ``GT=./., UV=0, SC=OctopuSV``).

Downstream VCF export may represent a true layout placeholder either as
missing (``./.``; the default) or operationally as reference (``0/0``) for
presence/absence cohort workflows.  The latter is an explicit export policy,
not a change to the SVCF fact layer.

Legacy sample blocks without SVCF 1.1 consensus fields are never reinterpreted.
Consensus itself remains implemented in ``sample_consensus``; this module only
interprets already-synthesized sample-mode blocks.
"""

from __future__ import annotations

from collections.abc import Mapping


UNOBSERVED_SAMPLE_GT_POLICIES = frozenset({"missing", "ref"})


def normalize_unobserved_sample_gt(value: object) -> str:
    """Validate and normalize the VCF export policy for unobserved samples."""
    policy = str(value).strip().lower()
    if policy not in UNOBSERVED_SAMPLE_GT_POLICIES:
        allowed = ", ".join(sorted(UNOBSERVED_SAMPLE_GT_POLICIES))
        raise ValueError(
            f"Invalid --unobserved-sample-gt value '{value}'. "
            f"Expected one of: {allowed}."
        )
    return policy


def _int_field(parsed_sample: Mapping[str, object], key: str):
    if key not in parsed_sample:
        return None
    raw = parsed_sample.get(key)
    if raw in (None, "", "."):
        return None
    try:
        return int(str(raw).strip())
    except (TypeError, ValueError):
        return None


def _is_missing_token(value: object) -> bool:
    return value in (None, "", ".", "unknown")


def is_unobserved_sample(parsed_sample: Mapping[str, object]) -> bool:
    """Return True only for the fixed-width *no-event* sample placeholder.

    ``UV=0`` alone is intentionally insufficient.  A real synthesized sample
    can have evidence yet contribute no valid vote (for example conflicting
    records from one caller); such a call must remain unresolved rather than
    being rewritten as reference under the ``ref`` export policy.

    The current SVCF 1.1 writer's placeholder contract is:
      * GT=0/0
      * UC=0 and UV=0
      * no synthesized source/event identity (SC and ID are missing)

    Legacy layouts without UC/UV therefore cannot be mistaken for a new
    placeholder.
    """
    if "UV" not in parsed_sample or "UC" not in parsed_sample:
        return False

    if _int_field(parsed_sample, "UV") != 0:
        return False
    if _int_field(parsed_sample, "UC") != 0:
        return False

    gt = parsed_sample.get("GT")
    if str(gt) != "0/0":
        return False

    if not _is_missing_token(parsed_sample.get("SC")):
        return False
    if not _is_missing_token(parsed_sample.get("ID")):
        return False

    return True


def downstream_sample_gt(
    parsed_sample: Mapping[str, object],
    *,
    default: str = "./.",
    unobserved_sample_gt: str = "missing",
) -> str:
    """Return the downstream genotype for one sample-mode SVCF block.

    Only a true SVCF 1.1 no-event placeholder is affected by
    ``unobserved_sample_gt``:

    * ``missing`` -> ``./.`` (default; does not assert homozygous reference)
    * ``ref``     -> ``0/0`` (explicit presence/absence cohort policy)

    Evidence-backed calls, unresolved calls with evidence, malformed fields,
    and legacy layouts are preserved exactly.
    """
    policy = normalize_unobserved_sample_gt(unobserved_sample_gt)

    raw_gt = parsed_sample.get("GT")
    gt = default if raw_gt in (None, "") else str(raw_gt)

    if not is_unobserved_sample(parsed_sample):
        return gt

    if policy == "ref":
        return "0/0"
    return "./."
