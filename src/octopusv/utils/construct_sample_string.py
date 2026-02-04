def construct_sample_string(sv_event):
    format_parts = sv_event.format.split(":")
    sample_parts = sv_event.sample.split(":")
    format_sample_dict = dict(zip(format_parts, sample_parts, strict=False))

    # Deal with AD, we have to calculate it by ourself.
    if "AD" in format_sample_dict:
        ad = format_sample_dict["AD"]
    elif "DR" in format_sample_dict and "DV" in format_sample_dict:
        dr = format_sample_dict["DR"]
        dv = format_sample_dict["DV"]
        ad = f"{dr},{dv}"
    else:
        ad = ".,."

    # Deal with others
    gt = format_sample_dict.get("GT", ".")

    svlen = sv_event.info.get("SVLEN", "0")  # Default to '0' if not present
    try:
        ln = str(abs(int(svlen)))  # Convert to absolute value string if possible
    except ValueError:
        ln = "."  # Use '.' if svlen is not an integer

    st = sv_event.info.get("STRAND", ".")
    qv = sv_event.qual
    ty = sv_event.info.get("SVTYPE", ".")
    id_ = sv_event.id
    sc = sv_event.source  # The source is a dynamic attribute added in parse_vcf()
    ref = sv_event.ref
    alt = sv_event.orig_alt
    chr2 = sv_event.info.get("CHR2", ".")
    end = sv_event.info.get("END", ".")

    # Calculate the correct end position for insertions
    if ty == "INS" and svlen.isdigit():
        end = sv_event.pos + int(svlen)  # Adjust end based on the length of the insertion

    # Build CO
    co = f"{sv_event.chrom}_{sv_event.pos}-{chr2}_{end}" if chr2 != "." and end != "." else "."

    return f"{gt}:{ad}:{ln}:{st}:{qv}:{ty}:{id_}:{sc}:{ref}:{alt}:{co}"
