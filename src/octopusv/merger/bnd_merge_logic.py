import re


def should_merge_bnd(event1, event2, delta=50):
    """Determines whether two BND events should be merged based on strict pattern matching.

    BND events require nearly identical patterns to be merged, as they represent
    unclassified breakpoints that should be conservatively handled.

    Args:
        event1: First BND event
        event2: Second BND event
        delta: Position uncertainty threshold (default: 50bp)

    Returns:
        bool: True if events should be merged, False otherwise
    """
    # Parse BND patterns and coordinates from both events
    pattern1, target_chr1, target_pos1 = parse_bnd_alt(event1.alt)
    pattern2, target_chr2, target_pos2 = parse_bnd_alt(event2.alt)

    # BND pattern must be identical (no reciprocal matching like TRA)
    if pattern1 != pattern2:
        return False

    # Source and target chromosomes must match exactly
    if event1.chrom != event2.chrom or target_chr1 != target_chr2:
        return False

    # Check position proximity within delta threshold
    start_diff = abs(event1.pos - event2.pos)
    target_diff = abs(target_pos1 - target_pos2)

    return start_diff <= delta and target_diff <= delta


def parse_bnd_alt(alt_field):
    """Parse BND ALT field to extract pattern, target chromosome, and position.

    Args:
        alt_field (str): BND ALT field value (e.g., "]chr2:186737138]N")

    Returns:
        tuple: (pattern, target_chromosome, target_position)
            - pattern: BND orientation pattern (e.g., "]p]t", "t[p[")
            - target_chromosome: Target chromosome name
            - target_position: Target position as integer
    """
    if not alt_field or alt_field == ".":
        return "UNKNOWN", None, None

    try:
        # Extract chromosome and position using regex
        # Pattern: ]chr2:186737138]N or N[chr2:186737138[
        match = re.search(r'([^\[\]]*)([\[\]])([^:\[\]]+):(\d+)([\[\]])([^\[\]]*)', alt_field)
        if not match:
            return "UNKNOWN", None, None

        prefix, bracket1, target_chr, target_pos_str, bracket2, suffix = match.groups()
        target_pos = int(target_pos_str)

        # Classify BND pattern based on bracket orientation and position
        if prefix and not suffix:  # Starts with sequence: N]chr:pos] or N[chr:pos[
            if bracket1 == ']' and bracket2 == ']':
                pattern = "t]p]"
            elif bracket1 == '[' and bracket2 == '[':
                pattern = "t[p["
            else:
                pattern = "MIXED"
        elif suffix and not prefix:  # Ends with sequence: ]chr:pos]N or [chr:pos[N
            if bracket1 == ']' and bracket2 == ']':
                pattern = "]p]t"
            elif bracket1 == '[' and bracket2 == '[':
                pattern = "[p[t"
            else:
                pattern = "MIXED"
        else:
            pattern = "UNKNOWN"

        return pattern, target_chr, target_pos

    except (ValueError, AttributeError, IndexError):
        return "UNKNOWN", None, None


def classify_bnd_pattern(alt_field):
    """Classify BND pattern from ALT field for simpler pattern comparison.

    Args:
        alt_field (str): BND ALT field value

    Returns:
        str: Simplified pattern classification
    """
    pattern, _, _ = parse_bnd_alt(alt_field)
    return pattern


def are_identical_bnd_patterns(pattern1, pattern2):
    """Check if two BND patterns are identical.

    Unlike TRA events, BND events must have identical patterns to be merged,
    as they represent unclassified breakpoints.

    Args:
        pattern1 (str): First BND pattern
        pattern2 (str): Second BND pattern

    Returns:
        bool: True if patterns are identical
    """
    return pattern1 == pattern2 and pattern1 != "UNKNOWN"