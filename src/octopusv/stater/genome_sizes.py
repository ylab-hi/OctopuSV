"""Built-in standard-chromosome lengths and contig-name handling for density.

Scope is deliberately limited to the standard chromosomes (1-22, X, Y, MT).
Decoy / unplaced / alt contigs (hs37d5, GL000*, NC_*, KI270*, ...) are NOT
given lengths: density on them is not meaningful, so they are reported as
unmatched rather than guessed. We never fall back to a length of 1.

Length sources, in priority order, are decided by the caller:
    --fai  (most trustworthy)  >  --genome  >  auto

GRCh37/hg19 and GRCh38/hg38 standard-chromosome lengths are public constants.
"""

# Standard chromosome lengths (bp). Keys are the bare names without 'chr'.
HG19_SIZES = {
    "1": 249250621, "2": 243199373, "3": 198022430, "4": 191154276,
    "5": 180915260, "6": 171115067, "7": 159138663, "8": 146364022,
    "9": 141213431, "10": 135534747, "11": 135006516, "12": 133851895,
    "13": 115169878, "14": 107349540, "15": 102531392, "16": 90354753,
    "17": 81195210, "18": 78077248, "19": 59128983, "20": 63025520,
    "21": 48129895, "22": 51304566, "X": 155270560, "Y": 59373566,
    "MT": 16569,
}

HG38_SIZES = {
    "1": 248956422, "2": 242193529, "3": 198295559, "4": 190214555,
    "5": 181538259, "6": 170805979, "7": 159345973, "8": 145138636,
    "9": 138394717, "10": 133797422, "11": 135086622, "12": 133275309,
    "13": 114364328, "14": 107043718, "15": 101991189, "16": 90338345,
    "17": 83257441, "18": 80373285, "19": 58617616, "20": 64444167,
    "21": 46709983, "22": 50818468, "X": 156040895, "Y": 57227415,
    "MT": 16569,
}

# Aliases for the same name in different styles -> the bare key above.
_GENOMES = {
    "hg19": HG19_SIZES, "grch37": HG19_SIZES,
    "hg38": HG38_SIZES, "grch38": HG38_SIZES,
}


def normalize_contig_name(name):
    """Map a contig name to a bare standard-chromosome key, or None.

    Handles only standard chromosomes:
        chr1 -> 1,  X/chrX -> X,  chrM/M/MT -> MT.
    Decoy / unplaced / alt contigs return None (they are not standard).
    """
    n = name[3:] if name.lower().startswith("chr") else name
    if n in ("M", "MT", "m"):
        return "MT"
    # Accept 1..22, X, Y only.
    if n in {str(i) for i in range(1, 23)} or n in ("X", "Y"):
        return n
    return None


def load_genome_sizes(genome):
    """Return the length dict for a named built-in genome (case-insensitive)."""
    key = genome.lower()
    if key not in _GENOMES:
        raise ValueError(
            f"Unknown genome '{genome}'. Available: {sorted(_GENOMES.keys())}."
        )
    return _GENOMES[key]


def load_fai_lengths(fai_path):
    """Load contig lengths from a .fai index file: {name: length}."""
    lengths = {}
    with open(fai_path) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 2:
                try:
                    lengths[parts[0]] = int(parts[1])
                except ValueError:
                    continue
    return lengths


def auto_detect_genome(contig_names):
    """Guess a built-in genome from the observed contig names.

    Returns (genome_name, assumed_bool):
      * 'hs37d5' present            -> ('hg19', False)  GRCh37/1000G reference
      * KI270*/chrUn_* present      -> ('hg38', False)  GRCh38-style alt contigs
      * only ambiguous primaries    -> ('hg38', True)   assumed, flagged
    """
    names = set(contig_names)
    if any(n == "hs37d5" or n.startswith("GL000") for n in names):
        return "hg19", False
    if any("KI270" in n or n.startswith("chrUn_") for n in names):
        return "hg38", False
    # Only 1/2/3.. or chr1/chr2.. with no decoy clue: cannot tell hg19 vs hg38.
    return "hg38", True


def resolve_lengths(contig_names, fai=None, genome="auto"):
    """Decide which lengths to use and describe the source.

    Priority: fai > explicit genome > auto-detect.
    Returns (lengths_dict, source_dict) where source_dict is:
        {"kind": "fai"|"builtin", "genome": <name or None>,
         "path": <fai path or None>, "assumed": bool}
    """
    if fai:
        return load_fai_lengths(fai), {
            "kind": "fai", "genome": None, "path": fai, "assumed": False,
        }
    if genome and genome.lower() != "auto":
        return load_genome_sizes(genome), {
            "kind": "builtin", "genome": genome.lower(), "path": None,
            "assumed": False,
        }
    # auto
    detected, assumed = auto_detect_genome(contig_names)
    return load_genome_sizes(detected), {
        "kind": "builtin", "genome": detected, "path": None, "assumed": assumed,
    }
