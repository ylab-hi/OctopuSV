"""Genome-wide SV Circos plotter for OctopuSV.

Draws a two-layer Circos overview from parsed SVCFEvent objects:
  - inner LINK layer: one arc per displayed SV, colored by SV type
  - outer DENSITY layer: histogram of the two breakpoints of each displayed SV

Both layers are built from the SAME display-filtered SV set, so the outer
density and the inner links always refer to the same events.

INS is excluded by default (single effective breakpoint, not a rearrangement
link); it can be folded into the density layer only via include_ins=True.
"""

import logging
import math
from collections import defaultdict

logger = logging.getLogger(__name__)

# hg38 primary chromosome sizes (fallback when no .fai / contig lengths given)
HG38_CHROM_SIZES = {
    "chr1": 248956422, "chr2": 242193529, "chr3": 198295559, "chr4": 190214555,
    "chr5": 181538259, "chr6": 170805979, "chr7": 159345973, "chr8": 145138636,
    "chr9": 138394717, "chr10": 133797422, "chr11": 135086622, "chr12": 133275309,
    "chr13": 114364328, "chr14": 107043718, "chr15": 101991189, "chr16": 90338345,
    "chr17": 83257441, "chr18": 80373285, "chr19": 58617616, "chr20": 64444167,
    "chr21": 46709983, "chr22": 50818468, "chrX": 156040895, "chrY": 57227415,
}
CHROM_ORDER = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]

SV_COLORS = {
    "DEL": "#2E6FB6",   # blue
    "DUP": "#2CA05A",   # green
    "INV": "#7A4FA3",   # purple
    "TRA": "#D6353D",   # red
}
INTRA_TYPES = ("DEL", "DUP", "INV")


def normalize_chrom(chrom):
    """Normalize a chromosome name to 'chrN' form; return as-is if unrecognized."""
    if chrom is None:
        return None
    chrom = str(chrom)
    if chrom.startswith("chr"):
        return chrom
    if chrom in [str(i) for i in range(1, 23)] + ["X", "Y", "M", "MT"]:
        return "chr" + chrom.replace("MT", "M")
    return chrom


def load_chrom_sizes_from_fai(fai_path):
    """Load primary chr1-22/X/Y sizes from a reference .fai index."""
    sizes = {}
    with open(fai_path) as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 2:
                continue
            c = normalize_chrom(parts[0])
            if c in CHROM_ORDER:
                sizes[c] = int(parts[1])
    return {c: sizes[c] for c in CHROM_ORDER if c in sizes}


class CircosPlotter:
    """Build a genome-wide SV Circos plot from a list of SVCFEvent objects.

    Args:
        events: list of SVCFEvent (from SVCFFileEventCreator).
        chrom_sizes: optional dict {chrom: length}. Defaults to hg38.
        require_pass: keep only FILTER == 'PASS' records (default False;
            OctopuSV merges callers whose FILTER is not standard PASS).
        tra_support: minimum SUPPORT for TRA (interchromosomal) links.
        intra_support: minimum SUPPORT for intra DEL/DUP/INV links.
        intra_min_span / intra_max_span: span window (bp) for intra links.
            Intra events larger than intra_max_span go to the oversized table.
        include_ins: if True, INS breakpoints are added to the DENSITY layer
            only (never drawn as links).
        ins_support: minimum SUPPORT for INS when include_ins is True.
        bin_size: density histogram bin size (bp).
        drop_del / drop_dup / drop_inv / no_tra: per-type on/off switches.
    """

    def __init__(
        self,
        events,
        chrom_sizes=None,
        require_pass=False,
        tra_support=3,
        intra_support=3,
        intra_min_span=5_000,
        intra_max_span=50_000_000,
        include_ins=False,
        ins_support=3,
        bin_size=5_000_000,
        drop_del=False,
        drop_dup=False,
        drop_inv=False,
        no_tra=False,
    ):
        self.events = events
        self.chrom_sizes = chrom_sizes or HG38_CHROM_SIZES.copy()
        self.require_pass = require_pass
        self.tra_support = tra_support
        self.intra_support = intra_support
        self.intra_min_span = intra_min_span
        self.intra_max_span = intra_max_span
        self.include_ins = include_ins
        self.ins_support = ins_support
        self.bin_size = bin_size
        self.drop_del = drop_del
        self.drop_dup = drop_dup
        self.drop_inv = drop_inv
        self.no_tra = no_tra

        self.links = []          # list of dicts: chrom1,pos1,chrom2,pos2,svtype,support
        self.breakpoints = []    # list of (chrom, pos) for density
        self.oversized = []      # intra events > intra_max_span
        self.stats = defaultdict(int)

    # ----------------------- extraction -----------------------
    @staticmethod
    def _support(event):
        val = event.info.get("SUPPORT")
        if val in (None, ".", ""):
            return None
        try:
            return int(float(val))
        except (TypeError, ValueError):
            return None

    @staticmethod
    def _span(event):
        svlen = event.info.get("SVLEN")
        if svlen not in (None, ".", ""):
            try:
                return abs(int(float(svlen)))
            except (TypeError, ValueError):
                pass
        if event.start_chrom == event.end_chrom:
            return abs(event.end_pos - event.start_pos)
        return None

    def _on_primary(self, chrom):
        return chrom in self.chrom_sizes

    def extract(self):
        """Populate self.links / self.breakpoints / self.oversized / self.stats."""
        for ev in self.events:
            self.stats["records"] += 1
            svtype = (ev.sv_type or "").upper()
            if svtype == "BND":
                svtype = "TRA"

            if self.require_pass and ev.filter != "PASS":
                self.stats["drop_not_pass"] += 1
                continue

            c1 = normalize_chrom(ev.start_chrom)
            c2 = normalize_chrom(ev.end_chrom)
            p1 = ev.start_pos
            p2 = ev.end_pos
            support = self._support(ev)

            # ---- INS: density-only, optional ----
            if svtype == "INS":
                if not self.include_ins:
                    self.stats["drop_ins"] += 1
                    continue
                if support is None or support < self.ins_support:
                    self.stats["drop_ins_lowsupport"] += 1
                    continue
                if self._on_primary(c1) and 1 <= p1 <= self.chrom_sizes[c1]:
                    self.breakpoints.append((c1, p1))
                    self.stats["ins_density"] += 1
                else:
                    self.stats["drop_ins_offchrom"] += 1
                continue

            # ---- type on/off switches ----
            if svtype == "TRA" and self.no_tra:
                self.stats["drop_tra_off"] += 1
                continue
            if svtype == "DEL" and self.drop_del:
                self.stats["drop_del_off"] += 1
                continue
            if svtype == "DUP" and self.drop_dup:
                self.stats["drop_dup_off"] += 1
                continue
            if svtype == "INV" and self.drop_inv:
                self.stats["drop_inv_off"] += 1
                continue

            is_inter = (svtype == "TRA") or (c1 != c2)

            # ---- TRA / interchromosomal ----
            if svtype == "TRA":
                if support is None or support < self.tra_support:
                    self.stats["drop_tra_lowsupport"] += 1
                    continue
                if not (self._on_primary(c1) and self._on_primary(c2)):
                    self.stats["drop_tra_offchrom"] += 1
                    continue
                if not (1 <= p1 <= self.chrom_sizes[c1] and 1 <= p2 <= self.chrom_sizes[c2]):
                    self.stats["drop_tra_badpos"] += 1
                    continue
                self.links.append(dict(svtype="TRA", chrom1=c1, pos1=p1,
                                       chrom2=c2, pos2=p2, support=support))
                self.breakpoints.append((c1, p1))
                self.breakpoints.append((c2, p2))
                self.stats["tra_link"] += 1
                continue

            # ---- intra DEL/DUP/INV ----
            if svtype in INTRA_TYPES:
                if is_inter or not self._on_primary(c1):
                    self.stats["drop_intra_offchrom"] += 1
                    continue
                span = self._span(ev)
                if span is None:
                    self.stats["drop_intra_nospan"] += 1
                    continue
                if span > self.intra_max_span:
                    self.oversized.append(dict(sv_id=ev.sv_id, svtype=svtype, chrom=c1,
                                               pos1=p1, pos2=p2, span=span,
                                               support=support, filter=ev.filter))
                    self.stats["intra_oversized"] += 1
                    continue
                if span < self.intra_min_span:
                    self.stats["drop_intra_small"] += 1
                    continue
                if support is None or support < self.intra_support:
                    self.stats["drop_intra_lowsupport"] += 1
                    continue
                if not (1 <= p1 <= self.chrom_sizes[c1] and 1 <= p2 <= self.chrom_sizes[c1]):
                    self.stats["drop_intra_badpos"] += 1
                    continue
                self.links.append(dict(svtype=svtype, chrom1=c1, pos1=p1,
                                       chrom2=c1, pos2=p2, support=support))
                self.breakpoints.append((c1, p1))
                self.breakpoints.append((c1, p2))
                self.stats[f"{svtype.lower()}_link"] += 1
                continue

            self.stats[f"drop_othertype_{svtype}"] += 1

    # ----------------------- density binning -----------------------
    def _bin_breakpoints(self):
        per_chrom = defaultdict(list)
        for c, p in self.breakpoints:
            per_chrom[c].append(p)
        bins = {}
        for chrom, size in self.chrom_sizes.items():
            n_bins = max(1, math.ceil(size / self.bin_size))
            counts = [0] * n_bins
            for p in per_chrom.get(chrom, []):
                idx = int((p - 1) // self.bin_size)
                if 0 <= idx < n_bins:
                    counts[idx] += 1
            bins[chrom] = counts
        return bins

    # ----------------------- plotting -----------------------
    def plot(
        self,
        output_file,
        sample_name="sample",
        tra_alpha=0.30,
        intra_alpha=0.65,
        link_lw=0.6,
        tra_lw=1.1,
        intra_height=0.30,
        tra_height=0.55,
        dpi=300,
    ):
        """Render the Circos figure to output_file (.pdf/.png/.svg)."""
        # Imports kept local so importing this module never hard-requires
        # matplotlib/pycirclize unless the user actually draws.
        import matplotlib.pyplot as plt
        from matplotlib.lines import Line2D
        from pycirclize import Circos

        if not self.links and not self.breakpoints:
            logger.warning("No SV passed the display filters; nothing to plot.")

        bins = self._bin_breakpoints()
        max_count = max((max(v) if v else 0 for v in bins.values()), default=0)
        max_count = max(1, max_count)

        circos = Circos(sectors=self.chrom_sizes, space=2)
        chr_r = (97, 100)
        dens_r = (84, 95)
        link_r = 82

        for sector in circos.sectors:
            chrom = sector.name
            ct = sector.add_track(chr_r)
            ct.axis(fc="#D9D9D9", ec="white", lw=0.5)
            sector.text(chrom.replace("chr", ""), r=104, size=8, orientation="vertical")
            if sector.size >= 60_000_000:
                ct.xticks_by_interval(50_000_000, tick_length=1.2, show_label=False,
                                      line_kws=dict(color="#888888", lw=0.4))
            dt = sector.add_track(dens_r)
            dt.axis(fc="white", ec="#E0E0E0", lw=0.3)
            counts = bins.get(chrom, [])
            if counts:
                size = sector.size
                mids = []
                for i in range(len(counts)):
                    start = i * self.bin_size
                    end = min((i + 1) * self.bin_size, size)
                    mids.append((start + end) / 2)
                dt.bar(mids, counts, width=self.bin_size * 0.9, vmin=0, vmax=max_count,
                       color="#4D4D4D", ec="#4D4D4D", lw=0, align="center")

        def draw(subset, alpha, lw, height):
            for r in subset:
                color = SV_COLORS.get(r["svtype"], "#999999")
                p1, p2 = int(r["pos1"]), int(r["pos2"])
                circos.link(
                    (r["chrom1"], max(1, p1 - 1), min(self.chrom_sizes[r["chrom1"]], p1 + 1)),
                    (r["chrom2"], max(1, p2 - 1), min(self.chrom_sizes[r["chrom2"]], p2 + 1)),
                    r1=link_r, r2=link_r, color=color, alpha=alpha, lw=lw, height_ratio=height,
                )

        # TRA first (bottom), intra on top
        draw([r for r in self.links if r["svtype"] == "TRA"], tra_alpha, tra_lw, tra_height)
        draw([r for r in self.links if r["svtype"] in INTRA_TYPES], intra_alpha, link_lw, intra_height)

        n = len(self.links)
        n_tra = sum(1 for r in self.links if r["svtype"] == "TRA")
        n_intra = n - n_tra
        n_bp = len(self.breakpoints)
        circos.text(
            f"{sample_name}\n{n:,} links\n(intra {n_intra:,} / TRA {n_tra:,})\n"
            f"{n_bp:,} breakpoints\nbin={self.bin_size / 1e6:g} Mb",
            r=22, size=9, ha="center", va="center",
        )

        fig = circos.plotfig()
        present = [t for t in ("DEL", "DUP", "INV", "TRA")
                   if any(r["svtype"] == t for r in self.links)]
        handles = [Line2D([0], [0], color=SV_COLORS[t], lw=2, label=t) for t in present]
        handles.append(Line2D([0], [0], color="#4D4D4D", lw=4,
                              label=f"{self.bin_size / 1e6:g} Mb breakpoint count"))
        circos.ax.legend(handles=handles, loc="lower center", bbox_to_anchor=(0.5, -0.07),
                         ncol=max(1, len(handles)), fontsize=8, frameon=False)
        fig.savefig(output_file, dpi=dpi, bbox_inches="tight")
        plt.close(fig)

    def write_oversized_table(self, path):
        """Write the intra > max-span events (inspect manually) to a TSV."""
        cols = ["sv_id", "svtype", "chrom", "pos1", "pos2", "span", "support", "filter"]
        with open(path, "w") as f:
            f.write("\t".join(cols) + "\n")
            for r in self.oversized:
                f.write("\t".join(str(r[c]) for c in cols) + "\n")