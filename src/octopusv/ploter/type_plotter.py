"""SV type distribution plotter.

Preferred input is the structured type stats dict from SVStater:

    {
        "counts": {"DEL": 26886, "DUP": 10563, ...},
        "percentages": {"DEL": 44.41, ...}
    }

For backward compatibility, this plotter also accepts full stat JSON or a
legacy stat.txt path.
"""

from __future__ import annotations

import collections
import json
import logging
from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt

LOGGER = logging.getLogger(__name__)

COLOR_SCHEME = {
    "TRA": "#FF6B6B",
    "INV": "#4ECDC4",
    "DUP": "#45B7D1",
    "INS": "#96CEB4",
    "DEL": "#6C5B7B",
}

SV_ORDER = ["TRA", "INV", "DUP", "INS", "DEL"]


def _read_json(path: Path) -> dict[str, Any]:
    with path.open() as handle:
        return json.load(handle)


class TypePlotter:
    """Plot SV type distribution."""

    def __init__(self, stats_or_file: dict[str, Any] | str | Path):
        self.stats_or_file = stats_or_file
        self.data = self._load_data(stats_or_file)

    def _load_data(self, stats_or_file: dict[str, Any] | str | Path) -> dict[str, tuple[int, float]]:
        if isinstance(stats_or_file, dict):
            return self._from_dict(stats_or_file)

        path = Path(stats_or_file)
        if path.suffix.lower() == ".json":
            return self._from_dict(_read_json(path))

        return self._from_legacy_text(path)

    def _from_dict(self, data: dict[str, Any]) -> dict[str, tuple[int, float]]:
        if "svtype_counts" in data:
            counts = data.get("svtype_counts", {}) or {}
            percentages = data.get("svtype_percentages", {}) or {}
        else:
            type_stats = data.get("type", data)
            counts = type_stats.get("counts", {}) or {}
            percentages = type_stats.get("percentages", {}) or {}

        total = sum(int(v) for v in counts.values())
        parsed: dict[str, tuple[int, float]] = {}

        for sv_type, count in counts.items():
            c = int(count)
            pct = percentages.get(sv_type)
            if pct is None:
                pct = (c / total * 100) if total else 0.0
            parsed[str(sv_type)] = (c, float(pct))

        return parsed

    def _from_legacy_text(self, path: Path) -> dict[str, tuple[int, float]]:
        sv_types: dict[str, tuple[int, float]] = {}
        in_section = False

        with path.open() as handle:
            for line in handle:
                if "SV Type Analysis" in line:
                    in_section = True
                    continue
                if in_section and not line.strip():
                    break
                if not in_section:
                    continue

                parts = line.strip().split("=")
                if len(parts) != 2:
                    continue

                sv_type = parts[0].strip()
                fields = parts[1].strip().split()
                if not fields:
                    continue

                try:
                    count = int(fields[0])
                    percentage = float(fields[1].strip("()%")) if len(fields) >= 2 else 0.0
                except ValueError:
                    continue

                sv_types[sv_type] = (count, percentage)

        return sv_types

    def _sort_data(self) -> collections.OrderedDict[str, tuple[int, float]]:
        ordered = collections.OrderedDict()

        for sv_type in SV_ORDER:
            if sv_type in self.data:
                ordered[sv_type] = self.data[sv_type]

        for sv_type, value in self.data.items():
            if sv_type not in ordered:
                ordered[sv_type] = value

        return ordered

    def plot(self, output_prefix: str | Path, *, save_svg: bool = True) -> None:
        """Create and save the SV type distribution plot."""
        if not self.data:
            LOGGER.error("No SV type data to plot.")
            return

        ordered_data = self._sort_data()
        types = list(ordered_data.keys())
        counts = [count for count, _ in ordered_data.values()]

        if sum(counts) == 0:
            LOGGER.error("All SV type counts are zero; nothing to plot.")
            return

        colors = [COLOR_SCHEME.get(sv_type, "#7f8c8d") for sv_type in types]

        fig, ax = plt.subplots(figsize=(10, 8))
        plt.rcParams["font.family"] = "sans-serif"
        plt.rcParams["svg.fonttype"] = "none"

        wedges, texts, autotexts = ax.pie(
            counts,
            labels=types,
            colors=colors,
            autopct="%1.1f%%",
            pctdistance=0.85,
            startangle=90,
            wedgeprops={"width": 0.5, "edgecolor": "white", "linewidth": 2},
        )

        plt.setp(autotexts, size=9, weight="bold", color="white")
        plt.setp(texts, size=10, weight="bold")

        centre_circle = plt.Circle((0, 0), 0.70, fc="white", ec="#e0e0e0")
        ax.add_artist(centre_circle)

        legend_labels = [
            f"{sv_type}: {ordered_data[sv_type][0]:,}"
            for sv_type in types
        ]
        ax.legend(
            wedges,
            legend_labels,
            title="SV Types",
            loc="center left",
            bbox_to_anchor=(1, 0, 0.5, 1),
            frameon=False,
            title_fontsize=12,
            fontsize=10,
        )

        ax.set_title("Structural Variant Type Distribution", pad=20, fontsize=14, fontweight="bold")
        ax.axis("equal")

        fig.tight_layout()
        output_prefix = str(output_prefix)
        fig.savefig(f"{output_prefix}.png", dpi=300, bbox_inches="tight", facecolor="white")

        if save_svg:
            fig.savefig(f"{output_prefix}.svg", bbox_inches="tight", facecolor="white")

        plt.close(fig)
        LOGGER.info("Type plot saved as %s.png%s", output_prefix, " and .svg" if save_svg else "")
