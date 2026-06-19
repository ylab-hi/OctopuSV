"""Stat analyzers that consume pre-parsed records and return structured dicts.

Every analyzer takes the shared `records` list (and `sample_names`) and returns
a plain dict of numbers -- no formatted strings, no regex round-trips. Text,
JSON and HTML are all rendered from these dicts downstream.

Correctness fixes baked in here:
  * SizeAnalyzer excludes TRA/BND (their END is a mate position, not a linear
    length) and records why each excluded record was dropped.
  * ChromosomeAnalyzer reports counts only; density requires real contig
    lengths (--fai) and is otherwise left out.
  * QCAnalyzer reports null (None) for missing QUAL/SUPPORT stats instead of
    a fake 0, so "no data" is distinguishable from "value is zero".
  * INFO parsing is centralized in stat_reader.parse_info.
"""

import statistics
from collections import Counter, defaultdict

from octopusv.utils.genotype_resolver import resolve_multi_caller_genotype


# ---------------------------------------------------------------------------
# Type
# ---------------------------------------------------------------------------

class TypeAnalyzer:
    """Count records per SVTYPE."""

    def __init__(self, records):
        self.records = records

    def analyze(self):
        counts = Counter(r.svtype for r in self.records)
        total = sum(counts.values())
        return {
            "total": total,
            "counts": dict(counts),
            "percentages": {
                t: (c / total * 100 if total else 0.0) for t, c in counts.items()
            },
        }


# ---------------------------------------------------------------------------
# Size
# ---------------------------------------------------------------------------

_NON_LINEAR = {"TRA", "BND"}

_SIZE_BINS = [
    ("0-50 bp", lambda s: s <= 50),
    ("51-100 bp", lambda s: s <= 100),
    ("101-500 bp", lambda s: s <= 500),
    ("501-1 kb", lambda s: s <= 1000),
    ("1 kb-10 kb", lambda s: s <= 10000),
    (">10 kb", lambda s: True),
]


class SizeAnalyzer:
    """Summaries of linear SV sizes, excluding TRA/BND."""

    def __init__(self, records, min_size=50, max_size=None):
        self.records = records
        self.min_size = min_size
        self.max_size = max_size

    def _linear_size(self, r):
        """Return an integer size for a linear SV, or None if not derivable."""
        svlen = r.info.get("SVLEN")
        if svlen not in (None, ".", True):
            try:
                return abs(int(svlen))
            except ValueError:
                pass
        # Fall back to END - POS for DEL/DUP/INV.
        end = r.info.get("END")
        try:
            return abs(int(end) - int(r.pos))
        except (ValueError, TypeError):
            return None

    def analyze(self):
        sizes = []
        bins = defaultdict(int)
        excluded = {
            "TRA_or_BND_no_linear_size": 0,
            "missing_or_invalid_size": 0,
            "below_min_size": 0,
            "above_max_size": 0,
        }

        for r in self.records:
            if r.svtype in _NON_LINEAR:
                excluded["TRA_or_BND_no_linear_size"] += 1
                continue
            size = self._linear_size(r)
            if size is None:
                excluded["missing_or_invalid_size"] += 1
                continue
            if size < self.min_size:
                excluded["below_min_size"] += 1
                continue
            if self.max_size is not None and size > self.max_size:
                excluded["above_max_size"] += 1
                continue
            sizes.append(size)
            for label, test in _SIZE_BINS:
                if test(size):
                    bins[label] += 1
                    break

        if sizes:
            summary = {
                "records_analyzed": len(sizes),
                "min": min(sizes),
                "max": max(sizes),
                "mean": statistics.mean(sizes),
                "median": statistics.median(sizes),
                "stdev": statistics.stdev(sizes) if len(sizes) > 1 else 0.0,
                "bins": dict(bins),
                "excluded": excluded,
            }
        else:
            summary = {
                "records_analyzed": 0,
                "min": None, "max": None, "mean": None,
                "median": None, "stdev": None,
                "bins": {}, "excluded": excluded,
            }
        return summary


# ---------------------------------------------------------------------------
# Chromosome
# ---------------------------------------------------------------------------

class ChromosomeAnalyzer:
    """Per-chromosome record counts, plus density when lengths are available.

    lengths: dict {bare_chrom_name: length_bp} for standard chromosomes only
             (from --fai or a built-in genome). May be None/empty.
    Density is computed only for contigs whose normalized name has a known
    length; everything else (decoy/unplaced/alt) goes to unmatched_contigs.
    There is no fallback to length 1.
    """

    def __init__(self, records, lengths=None):
        self.records = records
        self.lengths = lengths or {}

    def analyze(self):
        from octopusv.stater.genome_sizes import normalize_contig_name

        counts = Counter(r.chrom for r in self.records)
        result = {"counts": dict(counts)}

        if not self.lengths:
            result["density_per_mb"] = None
            result["unmatched_contigs"] = []
            return result

        density = {}
        unmatched = []
        for chrom, count in counts.items():
            norm = normalize_contig_name(chrom)
            length = self.lengths.get(norm) if norm else None
            if length:
                density[chrom] = count / (length / 1_000_000)
            else:
                unmatched.append(chrom)
        result["density_per_mb"] = density
        result["unmatched_contigs"] = unmatched
        return result


# ---------------------------------------------------------------------------
# QC
# ---------------------------------------------------------------------------

def _summary_or_none(values):
    """Quartile summary of a numeric list; all-None when the list is empty."""
    if not values:
        return {"n": 0, "mean": None, "median": None,
                "min": None, "max": None, "q1": None, "q3": None}
    q1 = statistics.quantiles(values, n=4)[0] if len(values) >= 4 else values[0]
    q3 = statistics.quantiles(values, n=4)[2] if len(values) >= 4 else values[-1]
    return {
        "n": len(values),
        "mean": statistics.mean(values),
        "median": statistics.median(values),
        "min": min(values),
        "max": max(values),
        "q1": q1,
        "q3": q3,
    }


class QCAnalyzer:
    """QUAL summary, FILTER distribution, SUPPORT summary."""

    def __init__(self, records):
        self.records = records

    def analyze(self):
        quals = []
        supports = []
        filters = Counter()

        for r in self.records:
            # QUAL: only count present numeric values.
            if r.qual not in (".", ""):
                try:
                    quals.append(float(r.qual))
                except ValueError:
                    pass
            filters[r.filter if r.filter else "UNKNOWN"] += 1
            support = r.info.get("SUPPORT")
            if support not in (None, ".", True):
                try:
                    supports.append(int(support))
                except ValueError:
                    pass

        total = sum(filters.values())
        return {
            "qual": _summary_or_none(quals),
            "support": _summary_or_none(supports),
            "filter_counts": dict(filters),
            "filter_percentages": {
                f: (c / total * 100 if total else 0.0) for f, c in filters.items()
            },
        }


# ---------------------------------------------------------------------------
# Genotype
# ---------------------------------------------------------------------------

class GenotypeAnalyzer:
    """Genotype distribution, mode-aware.

    Single/caller mode -> one resolved genotype per record (caller-mode uses
    the shared multi-caller voting rule). Sample/multi mode -> per-sample and
    overall distributions.
    """

    def __init__(self, records, sample_names):
        self.records = records
        self.sample_names = sample_names

    def _gt_from_segment(self, fmt, segment):
        keys = fmt.split(":")
        if "GT" not in keys:
            return None
        idx = keys.index("GT")
        parts = segment.split(":")
        return parts[idx] if idx < len(parts) else None

    def analyze(self):
        # Sample/multi mode: more than one declared sample column.
        if len(self.sample_names) > 1:
            return self._analyze_sample_mode()
        return self._analyze_caller_mode()

    def _analyze_caller_mode(self):
        genotypes = Counter()
        for r in self.records:
            if not r.sample_cols:
                continue
            if len(r.sample_cols) == 1:
                gt = self._gt_from_segment(r.format, r.sample_cols[0])
            else:
                # Reconstruct the INFO string for SOURCES-order tie-break.
                info_str = ";".join(
                    k if v is True else f"{k}={v}" for k, v in r.info.items()
                )
                gt = resolve_multi_caller_genotype(r.format, r.sample_cols, info_str)
            if gt:
                genotypes[gt] += 1
        total = sum(genotypes.values())
        return {
            "mode": "caller",
            "overall": dict(genotypes),
            "overall_percentages": {
                g: (c / total * 100 if total else 0.0) for g, c in genotypes.items()
            },
        }

    def _analyze_sample_mode(self):
        per_sample = {name: Counter() for name in self.sample_names}
        for r in self.records:
            cols = r.sample_cols[:len(self.sample_names)]
            for i, name in enumerate(self.sample_names):
                if i < len(cols):
                    gt = self._gt_from_segment(r.format, cols[i])
                    if gt:
                        per_sample[name][gt] += 1
        overall = Counter()
        for c in per_sample.values():
            overall.update(c)
        return {
            "mode": "sample",
            "per_sample": {name: dict(c) for name, c in per_sample.items()},
            "overall": dict(overall),
        }
