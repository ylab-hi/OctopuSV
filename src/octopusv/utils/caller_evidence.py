"""Caller-aware record support (INFO/SUPPORT) for ``octopusv correct``.

Each SV caller records its supporting evidence in its own field.  This module
holds one profile per caller whose field differs from the generic INFO keys
already understood by ``SVEvent`` (SUPPORT, SUPP, SUPPREAD, WEIGHT, RE).  The
support is resolved once, when the record is read, and every later step
(quality filtering, BND conversion, SVCF writing) uses that value.

Rules
-----
* A value is the caller's own evidence count for the alternative allele.
  Channels of one caller are added (e.g. Manta paired + split reads); this is
  an evidence count, not a count of distinct molecules.
* "Unknown" and "zero" are different: SUPPORT=0 only when the input states 0;
  a missing field, '.', or a negative sentinel (HitSV -1) is unknown.
* A channel missing in one sample does not make that sample unknown when the
  caller's other channel is present (Manta PR-only / SR-only records).
* Record support across samples: a record-level INFO total from the caller is
  used as is (LUMPY SU).  Otherwise per-sample values are summed only when every
  sample is known; if any sample is unknown the record support is unknown.
  Per-sample FORMAT fields are never modified here.
* Callers without a profile keep the historical behavior exactly.
"""

from __future__ import annotations

UNKNOWN = "."

# Profile name -> (human-readable source label, per-sample extractor or None,
# record-level INFO key or None).  Matching is in ``detect_caller``.
_SOURCE_LABELS = {
    "lumpy": "INFO/SU",
    "delly": "FORMAT/DV+RV",
    "manta": "FORMAT/PR_ALT+SR_ALT",
    "pbsv": "FORMAT/AD_ALT",
    "hitsv": "FORMAT/LR_ALT+SR_ALT",
    "svaba": "FORMAT/AD",
}


def detect_caller(source_info, info):
    """Return the profile name for a record, or None for the historical path.

    ``source_info`` is the ``##source`` header value; DELLY writes no
    ``##source`` line, so it is recognized by its per-record SVMETHOD.
    """
    source = (source_info or "").lower()
    method = str(info.get("SVMETHOD") or "").upper()
    if method.startswith("EMBL.DELLY") or source.startswith("delly"):
        return "delly"
    if source.startswith("lumpy"):
        return "lumpy"
    if source.startswith("generatesvcandidates") or source.startswith("manta"):
        return "manta"
    if source.startswith("pbsv"):
        return "pbsv"
    if source.startswith("hitsv"):
        return "hitsv"
    if source.startswith("svaba"):
        return "svaba"
    return None


def source_label(caller):
    return _SOURCE_LABELS.get(caller, "historical INFO keys / FORMAT PR")


def _count(value):
    """Non-negative integer, or None for missing, '.', negative sentinel or junk."""
    if value is None:
        return None
    text = str(value).strip()
    if text in ("", "."):
        return None
    try:
        number = int(text)
    except ValueError:
        try:
            number = float(text)
        except ValueError:
            return None
        if number != int(number):
            return None
        number = int(number)
    return number if number >= 0 else None


def _alt_of_ref_alt_list(value):
    """ALT total of a 'ref,alt[,alt...]' list (Number=R), or None."""
    if value is None:
        return None
    parts = str(value).split(",")
    if len(parts) < 2:
        return None
    alts = [_count(part) for part in parts[1:]]
    if any(item is None for item in alts):
        return None
    return sum(alts)


def _channels(*values):
    """Sum of the known channels; unknown only if every channel is unknown."""
    known = [value for value in values if value is not None]
    return sum(known) if known else None


def _sample_dict(format_field, raw_sample):
    keys = str(format_field or "").split(":")
    values = str(raw_sample or "").split(":")
    return dict(zip(keys, values, strict=False))


def _sample_support(caller, fields):
    if caller == "delly":
        return _channels(_count(fields.get("DV")), _count(fields.get("RV")))
    if caller == "manta":
        return _channels(_alt_of_ref_alt_list(fields.get("PR")),
                         _alt_of_ref_alt_list(fields.get("SR")))
    if caller == "pbsv":
        return _alt_of_ref_alt_list(fields.get("AD"))
    if caller == "hitsv":
        # LR / SR are "ref,alt,unknown"; -1 marks a channel without data.
        def alt(value):
            parts = str(value).split(",") if value is not None else []
            return _count(parts[1]) if len(parts) >= 2 else None
        return _channels(alt(fields.get("LR")), alt(fields.get("SR")))
    if caller == "svaba":
        # SvABA AD is Number=1: the reads supporting the variant.
        value = fields.get("AD")
        if value is None or "," in str(value):
            return None
        return _count(value)
    if caller == "lumpy":
        return _count(fields.get("SU"))
    return None


def sample_supports(caller, format_field, samples):
    """Per-sample support (int or None), in sample order."""
    return [_sample_support(caller, _sample_dict(format_field, raw)) for raw in samples]


def resolve_support(caller, info, format_field, samples):
    """Record-level support string for a profiled caller ('.' if unknown)."""
    return _resolve(caller, info, format_field, samples)[0]


def _resolve(caller, info, format_field, samples):
    """(record support string, per-sample list or None if from an INFO total)."""
    if caller == "lumpy":
        total = _count(info.get("SU"))
        if total is not None:
            return str(total), None
    per_sample = sample_supports(caller, format_field, samples)
    if not per_sample or any(value is None for value in per_sample):
        return UNKNOWN, per_sample
    return str(sum(per_sample)), per_sample


def apply_caller_support(event, source_info):
    """Resolve INFO/SUPPORT once on a freshly parsed SVEvent of a profiled caller.

    Records from callers without a profile are left untouched, so their
    historical SUPPORT parsing (SVEvent INFO keys, first-sample PR) and the
    historical QualityFilter fallbacks stay exactly as before.
    """
    caller = detect_caller(source_info, event.info)
    event.support_caller = caller
    if caller is None:
        return None
    event.info["SUPPORT"], event.support_per_sample = _resolve(
        caller, event.info, event.format, event.samples)
    event.support_resolved = True
    return caller
