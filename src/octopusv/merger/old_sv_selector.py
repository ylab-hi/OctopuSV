def get_caller_name(sample_dict):
    """Get standardized caller name from sample dictionary.
    Handles different formats of caller names.
    """
    caller = sample_dict.get("SC", "unknown")
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
    all_source_files = set()
    merged_samples_list = []

    for sv in sv_group:
        # Add source file
        source_files = sv.source_file.split(",")
        all_source_files.update(source_files)

        # Add sample information
        sample_info = (sv.sample_name, sv.format, sv.sample)
        merged_samples_list.append(sample_info)

        # Calculate support score
        support_value = sv.info.get("SUPPORT", "0")
        if support_value in [None, ".", ""]:
            support = 0
        else:
            try:
                support = int(support_value)
            except ValueError:
                support = 0
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
        if total_score > max_score:
            max_score = total_score
            representative_sv = sv
        elif total_score == max_score:
            # Tiebreaker 1: Higher read support
            curr_support_value = representative_sv.info.get("SUPPORT", "0")
            current_support = 0 if curr_support_value in [None, ".", ""] else int(curr_support_value)
            if support > current_support:
                representative_sv = sv
            elif support == current_support:
                # Tiebreaker 2: Higher quality score
                curr_qual_value = representative_sv.quality
                current_qual = 50 if curr_qual_value in [None, ".", ""] else float(curr_qual_value)
                if qual > current_qual:
                    representative_sv = sv

    if representative_sv:
        representative_sv.source_file = ",".join(sorted(all_source_files))
        representative_sv.merged_samples = merged_samples_list

    return representative_sv
