import logging


def get_caller_name(sample_dict):
    """Get standardized caller name from sample dictionary.
    Handles different formats of caller names.
    """
    caller = sample_dict.get("SC", "unknown")
    # DELLY writes no ##source; `correct` names it from INFO/SVMETHOD
    # (EMBL.DELLY...).  Its evidence used to carry SC=".", so DELLY was never
    # ranked in representative scoring.  It stays unranked: ranking it would
    # change representatives, and the scoring needs its own review first
    # (QUAL is not comparable between callers).
    if str(caller).upper().startswith("EMBL.DELLY"):
        return caller
    # Standardize caller names
    if "svim" in caller.lower():
        return "SVIM"
    if "pbsv" in caller.lower():
        return "PBSV"
    if "cutesv" in caller.lower():
        return "cuteSV"
    if "sniffles" in caller.lower():
        return "sniffles"
    if "debreak" in caller.lower():
        return "debreak"
    if "manta" in caller.lower():
        return "manta"
    if "delly" in caller.lower():
        return "delly"
    if "svaba" in caller.lower():
        return "svaba"
    if "svdss" in caller.lower():
        return "svdss"
    return caller


_warned_support_sources = set()


def _support_count(sv):
    """SUPPORT as an int for scoring; missing, "." or invalid scores as 0.

    SVCF 1.1 requires SUPPORT to be "." or a non-negative integer.  An invalid
    value (12.5, NA, -3, a bare flag) is a non-blocking validator error
    (E_SUPPORT_001), so such files still merge: the value is treated as
    unknown, exactly like ".", never reinterpreted, and one warning per
    source names the first offending record.
    """
    value = sv.info.get("SUPPORT", "0")
    if value in (None, ".", ""):
        return 0
    text = "" if value is True else str(value).strip()
    if text.isascii() and text.isdigit():
        return int(text)
    source = getattr(sv, "_original_source_file", None) or getattr(sv, "source_file", "?")
    if source not in _warned_support_sources:
        _warned_support_sources.add(source)
        logging.warning(
            "Invalid SUPPORT=%r (record %r, %s): SVCF 1.1 expects '.' or a non-negative "
            "integer; treated as unknown for representative selection. Further invalid "
            "SUPPORT values from this input are not reported.",
            value, getattr(sv, "sv_id", "?"), source,
        )
    return 0


def get_caller_score(caller_name, caller_priority_list):
    """Calculate a score for a structural variant caller based on its priority ranking."""
    try:
        rank = caller_priority_list.index(caller_name) + 1
    except ValueError:
        # If caller not in priority list, rank it last
        rank = len(caller_priority_list) + 1
    total_callers = len(caller_priority_list) + 1
    return (total_callers - rank + 1) / total_callers * 10


def select_representative_sv(sv_group, weights=None, max_support=30, max_qual=100):
    """Select a representative SV record from a group of overlapping structural variants
    using a comprehensive scoring system.
    """
    if weights is None:
        weights = {"support": 0.5, "quality": 0.3, "caller": 0.2}

    # Predefined priority list for SV callers
    caller_priority_list = [
        "SVIM",  # SVIM has highest priority
        "PBSV",
        "cuteSV",
        "sniffles",
        "debreak",
        "manta",
        "delly",
        "svaba",
        "svdss",
    ]

    max_score = -1
    representative_sv = None
    representative_support = 0
    representative_qual = 50.0
    all_source_files = set()
    merged_samples_list = []
    merged_sample_records = []

    for sv in sv_group:
        # Preserve the exact source attached to this individual evidence record.
        # select_representative_sv() mutates representative_sv.source_file into a
        # comma-separated merged-source string, and the same group can be selected
        # more than once during one merge run. Keep an immutable per-event copy so
        # later calls never lose the original provenance.
        original_source_file = getattr(sv, "_original_source_file", None)
        if original_source_file is None:
            original_source_file = str(sv.source_file)
            sv._original_source_file = original_source_file

        source_files = [
            source.strip()
            for source in original_source_file.split(",")
            if source.strip()
        ]
        all_source_files.update(source_files)

        # Keep the existing three-field structure unchanged for caller-mode and
        # backward compatibility with existing consumers.
        sample_info = (sv.sample_name, sv.format, sv.sample)
        merged_samples_list.append(sample_info)

        # Sample mode needs an explicit source -> evidence binding. Store that in
        # a separate four-field structure instead of guessing the source later
        # from SAMPLE, ID, SC, or event order.
        merged_sample_records.append(
            (original_source_file, sv.sample_name, sv.format, sv.sample)
        )

        # Calculate support score
        support = _support_count(sv)
        support_score = min(support / max_support, 1.0) * 10

        # Calculate quality score
        qual_value = sv.quality
        if qual_value in [None, ".", ""]:
            qual = 50  # Default quality score
        else:
            try:
                qual = float(qual_value)
            except ValueError:
                qual = 50
        quality_score = min(qual / max_qual, 1.0) * 10

        # Calculate caller score using standardized caller name
        caller_name = get_caller_name(sv.sample)
        caller_score = get_caller_score(caller_name, caller_priority_list)

        # Calculate weighted total score
        total_score = (
            weights["support"] * support_score + weights["quality"] * quality_score + weights["caller"] * caller_score
        )

        # Update representative SV if current SV has higher score
        # Tiebreakers reuse the values already scored for the current
        # representative instead of converting its fields a second time.
        if total_score > max_score:
            take = True
        elif total_score == max_score:
            # Tiebreaker 1: higher read support; 2: higher quality
            take = support > representative_support or (
                support == representative_support and qual > representative_qual
            )
        else:
            take = False
        if take:
            max_score = total_score
            representative_sv = sv
            representative_support = support
            representative_qual = qual

    if representative_sv:
        representative_sv.source_file = ",".join(sorted(all_source_files))
        representative_sv.merged_samples = merged_samples_list
        representative_sv.merged_sample_records = merged_sample_records

    return representative_sv
