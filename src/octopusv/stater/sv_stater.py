"""SVStater: run all analyzers once and render text / JSON / HTML from dicts.

Single source of truth is `self.stats` (nested dicts of numbers). Three
renderers consume it:
  * to_text()      -> the human report (replaces the old write_results()).
  * to_json_dict() -> structured result for agents (stat --json).
  * to_html_dict() -> adapts stats into the exact keys ReportGenerator's
                      Jinja template expects, so --report keeps working
                      without touching the template.

No analyzer returns strings and nothing is parsed back out of text.
"""

import json
from pathlib import Path

from .stat_reader import read_records
from .genome_sizes import resolve_lengths
from .stat_analyzers import (
    TypeAnalyzer, SizeAnalyzer, ChromosomeAnalyzer, QCAnalyzer, GenotypeAnalyzer,
)


class SVStater:
    def __init__(self, input_file, min_size=50, max_size=None, fai=None, genome="auto"):
        self.input_file = input_file
        self.min_size = min_size
        self.max_size = max_size
        self.fai = fai
        self.genome = genome
        self.records = []
        self.sample_names = []
        self.stats = {}
        self.length_source = None

    def analyze(self):
        """Read the file once and run every analyzer into self.stats."""
        self.records, self.sample_names = read_records(self.input_file)

        # Decide chromosome lengths: --fai > --genome > auto-detect.
        contigs = {r.chrom for r in self.records}
        lengths, self.length_source = resolve_lengths(
            contigs, fai=self.fai, genome=self.genome
        )

        self.stats = {
            "type": TypeAnalyzer(self.records).analyze(),
            "size": SizeAnalyzer(self.records, self.min_size, self.max_size).analyze(),
            "chromosome": ChromosomeAnalyzer(self.records, lengths).analyze(),
            "qc": QCAnalyzer(self.records).analyze(),
            "genotype": GenotypeAnalyzer(self.records, self.sample_names).analyze(),
        }

    # -- JSON ----------------------------------------------------------------

    def to_json_dict(self):
        """Agent-facing structured result."""
        t = self.stats["type"]
        s = self.stats["size"]
        c = self.stats["chromosome"]
        qc = self.stats["qc"]
        g = self.stats["genotype"]

        warnings = []
        if self.length_source and self.length_source.get("assumed"):
            warnings.append(
                f"reference_assumed_{self.length_source.get('genome')}_pass_--genome_or_--fai_to_be_explicit"
            )
        if c.get("unmatched_contigs"):
            warnings.append("density_not_computed_for_some_contigs_see_unmatched_contigs")

        return {
            "input": self.input_file,
            "parameters": {"min_size": self.min_size, "max_size": self.max_size},
            "records": {
                "total": len(self.records),
                "size_analyzed": s["records_analyzed"],
            },
            "svtype_counts": t["counts"],
            "size": {
                "min": s["min"], "max": s["max"], "mean": s["mean"],
                "median": s["median"], "stdev": s["stdev"],
                "bins": s["bins"], "excluded": s["excluded"],
            },
            "chromosome": {
                "counts": c["counts"],
                "density_per_mb": c["density_per_mb"],
                "length_source": self.length_source,
                "unmatched_contigs": c.get("unmatched_contigs", []),
            },
            "filter_counts": qc["filter_counts"],
            "quality": qc["qual"],
            "support": qc["support"],
            "genotypes": g,
            "warnings": warnings,
        }

    def write_json(self, output_file=None):
        payload = json.dumps(self.to_json_dict(), indent=2)
        if output_file is None:
            return payload
        with Path(output_file).open("w") as f:
            f.write(payload)
        return None

    # -- TEXT ----------------------------------------------------------------

    def to_text(self):
        """Human report. Format mirrors the old report closely so existing
        plotters that read stat.txt keep working; the only intentional
        differences are the removed fake density and TRA-excluded sizes."""
        t = self.stats["type"]
        s = self.stats["size"]
        c = self.stats["chromosome"]
        qc = self.stats["qc"]
        g = self.stats["genotype"]
        L = []

        L.append("OctopusV report")
        L.append("-" * 40)
        L.append("")
        L.append(">>>>>>> Input")
        L.append("")
        L.append(f"input file = {self.input_file}")
        L.append("")

        # Type
        L.append(">>>>>>> Type Analysis")
        L.append("")
        L.append("SV Type Analysis               = ")
        for sv_type, count in t["counts"].items():
            pct = t["percentages"].get(sv_type, 0.0)
            L.append(f"{sv_type:30} = {count} ({pct:.2f}%)")
        L.append("")

        # Size
        L.append(">>>>>>> Size Analysis")
        L.append("")
        L.append(f"{'Total SVs analyzed':30} = {s['records_analyzed']}")
        if s["records_analyzed"]:
            L.append(f"{'Minimum size':30} = {s['min']} bp")
            L.append(f"{'Maximum size':30} = {s['max']} bp")
            L.append(f"{'Mean size':30} = {s['mean']:.2f} bp")
            L.append(f"{'Median size':30} = {s['median']} bp")
            L.append(f"{'Standard deviation':30} = {s['stdev']:.2f} bp")
            L.append("")
            L.append("Size distribution              = ")
            for label, count in s["bins"].items():
                L.append(f"{label:30} = {count}")
        L.append("")

        # Chromosome
        L.append(">>>>>>> Chromosome Analysis")
        L.append("")
        src = self.length_source or {}
        if src.get("kind") == "fai":
            src_label = f"length_source=fai:{src.get('path')}"
        elif src.get("kind") == "builtin":
            src_label = f"length_source=builtin_{src.get('genome')}"
            if src.get("assumed"):
                src_label += ",assumed"
                L.append(f"# reference assumed {src.get('genome')}; "
                         f"pass --genome or --fai to be explicit")
        else:
            src_label = "length_source=none"
        L.append("Chromosome Distribution        = ")
        density = c.get("density_per_mb") or {}
        for chrom, count in c["counts"].items():
            if chrom in density:
                L.append(f"{chrom:30} = {count} SVs ({density[chrom]:.2f} SVs/Mb; {src_label})")
            else:
                L.append(f"{chrom:30} = {count} SVs (density not available: contig length not found)")
        L.append("")

        # QC
        L.append(">>>>>>> Qc Analysis")
        L.append("")
        L.append("Quality Score (QUAL) Statistics = ")
        q = qc["qual"]
        L.append(f"{'Number of variants with QUAL':30} = {q['n']}")
        L.append(f"{'Average QUAL':30} = {_fmt(q['mean'])}")
        L.append(f"{'Median QUAL':30} = {_fmt(q['median'])}")
        L.append(f"{'Min QUAL':30} = {_fmt(q['min'])}")
        L.append(f"{'Max QUAL':30} = {_fmt(q['max'])}")
        L.append("")
        L.append("Filter Status                  = ")
        total_f = sum(qc["filter_counts"].values())
        for status, count in sorted(qc["filter_counts"].items()):
            pct = (count / total_f * 100) if total_f else 0.0
            L.append(f"{status:30} = {count} ({pct:.2f}%)")
        L.append("")
        L.append("Read Support Statistics        = ")
        sup = qc["support"]
        L.append(f"{'Number of variants with support info':30} = {sup['n']}")
        L.append(f"{'Average Read Support':30} = {_fmt(sup['mean'])}")
        L.append(f"{'Median Read Support':30} = {_fmt(sup['median'])}")
        L.append("")

        # Genotype
        L.append(">>>>>>> Genotype Analysis")
        L.append("")
        L.append("Genotype Distribution          = ")
        if g.get("mode") == "sample":
            for name, dist in g["per_sample"].items():
                tot = sum(dist.values())
                L.append(f"{name} Genotypes:")
                for gt, count in dist.items():
                    pct = (count / tot * 100) if tot else 0.0
                    L.append(f"  {gt:28} = {count} ({pct:.2f}%)")
            tot = sum(g["overall"].values())
            L.append("Overall:")
            for gt, count in g["overall"].items():
                pct = (count / tot * 100) if tot else 0.0
                L.append(f"  {gt:28} = {count} ({pct:.2f}%)")
        else:
            tot = sum(g["overall"].values())
            for gt, count in g["overall"].items():
                pct = (count / tot * 100) if tot else 0.0
                L.append(f"{gt:30} = {count} ({pct:.2f}%)")
        L.append("")
        return "\n".join(L)

    def write_results(self, output_file):
        with Path(output_file).open("w") as f:
            f.write(self.to_text())

    # -- HTML adapter --------------------------------------------------------

    def to_html_dict(self):
        """Adapt stats into the keys ReportGenerator's template expects.

        The legacy export_html produced these keys; we reproduce them directly
        from structured data instead of parsing formatted text.
        """
        t = self.stats["type"]
        s = self.stats["size"]
        qc = self.stats["qc"]
        g = self.stats["genotype"]

        sv_types = {
            sv_type: (count, t["percentages"].get(sv_type, 0.0))
            for sv_type, count in t["counts"].items()
        }
        filter_status = {
            status: (str(count), qc["filter_percentages"].get(status, 0.0))
            for status, count in qc["filter_counts"].items()
        }
        # genotype_dist: caller mode -> overall; sample mode -> overall too.
        genotype_dist = {
            gt: (count, 0.0) for gt, count in g.get("overall", {}).items()
        }

        return {
            "input_file": self.input_file,
            "output_file": self.input_file,
            "sv_types": sv_types,
            "total_svs": s["records_analyzed"],
            "min_size": s["min"] or 0,
            "max_size": s["max"] or 0,
            "mean_size": s["mean"] or 0.0,
            "median_size": s["median"] or 0.0,
            "std_dev": s["stdev"] or 0.0,
            "size_distribution": s["bins"],
            "avg_qual": qc["qual"]["mean"] or 0.0,
            "filter_status": filter_status,
            "avg_read_support": qc["support"]["mean"] or 0.0,
            "genotype_dist": genotype_dist,
            "sample_genotypes": g.get("per_sample", {}),
            "population_genotypes": {},
        }

    # Back-compat alias: old code called export_html(output_file).
    def export_html(self, output_file=None):
        return self.to_html_dict()


def _fmt(v):
    """Format a number or 'NA' for None (missing), for the text report."""
    return "NA" if v is None else f"{v:.2f}"
