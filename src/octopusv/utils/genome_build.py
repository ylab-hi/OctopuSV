"""Genome build inference from declared ##contig lengths.

Used by `octopusv header` to add an `inferred_build` field to the declared
contract. Inference is based on chromosome LENGTHS (the reliable fingerprint),
never on names, and is deliberately conservative:

  * A build is reported only if every key contig present in the file agrees
    with that build's expected length.
  * Any mismatch -> build "unknown" with the conflicting contigs listed, so we
    never confidently mislabel a mixed or custom reference. Unknown is safer
    than a wrong call.

This module is self-contained: it carries its own minimal contig-id
normalization (chr-prefix only) so the header reader does not need to depend on
the filtering module.
"""

# Key contigs chosen for high discrimination across builds.
# Keyed by normalized contig id (no 'chr' prefix, 'M' for mitochondrion).
BUILD_SIGNATURES = {
    "GRCh37": {"1": 249250621, "2": 243199373, "X": 155270560},
    "GRCh38": {"1": 248956422, "2": 242193529, "X": 156040895},
    "CHM13v2": {"1": 248387328, "2": 242696752, "X": 154259566},
}


def normalize_contig_id(contig):
    """Minimal contig-id normalization for build lookup.

    Strips a leading 'chr' and maps mitochondrial aliases to 'M'. This is only
    for matching against BUILD_SIGNATURES keys; it does not touch file content.
    """
    if contig is None:
        return None
    cid = str(contig).strip()
    if cid.lower().startswith("chr"):
        cid = cid[3:]
    if cid in ("MT", "M", "mt", "m"):
        return "M"
    # Normalize X/Y casing.
    if cid.upper() in ("X", "Y"):
        return cid.upper()
    return cid


def infer_build(contig_lengths):
    """Infer genome build from a {contig_id: length} mapping.

    Args:
        contig_lengths: dict of raw contig id -> int length (from ##contig).

    Returns:
        dict with keys: build, confidence, evidence, conflicts.
            build       -- "GRCh37" / "GRCh38" / "CHM13v2" / "unknown"
            confidence  -- "high" when a build matched, else "none"
            evidence    -- {contig: length} actually used to decide
            conflicts   -- contigs whose length matched no/ wrong build
                           (only populated when build is "unknown")
    """
    # Normalize the incoming lengths once.
    norm_lengths = {}
    for cid, length in contig_lengths.items():
        ncid = normalize_contig_id(cid)
        if ncid is not None and length is not None:
            norm_lengths[ncid] = length

    # Try each candidate build: it wins only if every key contig we HAVE agrees.
    for build, signature in BUILD_SIGNATURES.items():
        checked = {}
        agree = True
        for cid, expected in signature.items():
            if cid in norm_lengths:
                checked[cid] = norm_lengths[cid]
                if norm_lengths[cid] != expected:
                    agree = False
        if checked and agree:
            return {
                "build": build,
                "confidence": "high",
                "evidence": checked,
                "conflicts": [],
            }

    # No build agreed. Report which key contigs we saw and what length they had,
    # so an agent can see why we couldn't decide.
    conflicts = []
    seen_evidence = {}
    key_contigs = set()
    for signature in BUILD_SIGNATURES.values():
        key_contigs.update(signature.keys())
    for cid in sorted(key_contigs):
        if cid in norm_lengths:
            seen_evidence[cid] = norm_lengths[cid]
            conflicts.append(cid)

    return {
        "build": "unknown",
        "confidence": "none",
        "evidence": seen_evidence,
        "conflicts": conflicts,
    }