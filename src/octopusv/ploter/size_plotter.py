# size_plotter.py

"""SV size distribution plotter.

Preferred input is the structured size stats dict from SVStater:

    {
        "bins": {
            "0-50 bp": 0,
            "51-100 bp": 8935,
            ...
        },
        ...
    }

For backward compatibility, this plotter also accepts full stat JSON or a
legacy stat.txt path.
"""

from __future__ import annotations

import json
import logging
from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt

LOGGER = logging.getLogger(__name__)

SIZE_ORDER = ["0-50 bp", "51-100 bp", "101-500 bp", "501-1 kb", "1 kb-10 kb", ">10 kb"]


def _read_json(path: Path) -> dict[str, Any]:
    with path.open() as handle:
        return json.load(handle)


class SizePlotter:
    """Plot SV size-bin distribution."""

    def __init__(self, stats_or_file: dict[str, Any] | str | Path):
        self.stats_or_file = stats_or_file
        self.data = self._load_data(stats_or_file)

    def _load_data(self, stats_or_file: dict[str, Any] | str | Path) -> dict[str, int]:
        if isinstance(stats_or_file, dict):
            return self._from_dict(stats_or_file)

        path = Path(stats_or_file)
        if path.suffix.lower() == ".json":
            return self._from_dict(_read_json(path))

        return self._from_legacy_text(path)

    def _from_dict(self, data: dict[str, Any]) -> dict[str, int]:
        size_stats = data.get("size", data)
        bins = size_stats.get("bins", {}) or {}
        return {str(k): int(v) for k, v in bins.items()}

    def _from_legacy_text(self, path: Path) -> dict[str, int]:
        bins: dict[str, int] = {}
        in_section = False

        with path.open() as handle:
            for line in handle:
                if "Size distribution" in line:
                    in_section = True
                    continue
                if in_section and not line.strip():
                    break
                if not in_section:
                    continue

                parts = line.strip().split("=")
                if len(parts) != 2:
                    continue

                label = parts[0].strip()
                try:
                    bins[label] = int(parts[1].strip())
                except ValueError:
                    continue

        return bins

    def plot(self, output_prefix: str | Path, *, save_svg: bool = True) -> None:
        """Create and save the SV size distribution plot."""
        if not self.data:
            LOGGER.error("No size distribution data to plot.")
            return

        sizes = [self.data.get(label, 0) for label in SIZE_ORDER]
        x = range(len(SIZE_ORDER))

        fig, ax = plt.subplots(figsize=(12, 7))
        bars = ax.bar(x, sizes, width=0.7, edgecolor="white", linewidth=1.5)

        ax.set_xlabel("SV Size Range", fontsize=12, labelpad=10, fontweight="bold")
        ax.set_ylabel("Count", fontsize=12, labelpad=10, fontweight="bold")
        ax.set_title("Structural Variant Size Distribution", fontsize=14, pad=20, fontweight="bold")
        ax.set_xticks(list(x))
        ax.set_xticklabels(SIZE_ORDER, rotation=25, ha="right")
        ax.grid(True, axis="y", linestyle="--", alpha=0.3)
        ax.set_axisbelow(True)

        for spine in ax.spines.values():
            spine.set_visible(True)
            spine.set_linewidth(0.8)

        for bar in bars:
            height = bar.get_height()
            if height == 0:
                continue
            ax.text(
                bar.get_x() + bar.get_width() / 2.0,
                height,
                f"{int(height):,}",
                ha="center",
                va="bottom",
                fontsize=10,
                fontweight="bold",
            )

        fig.tight_layout()
        output_prefix = str(output_prefix)
        fig.savefig(f"{output_prefix}.png", dpi=300, bbox_inches="tight", facecolor="white")

        if save_svg:
            fig.savefig(f"{output_prefix}.svg", format="svg", bbox_inches="tight", facecolor="white")

        plt.close(fig)
        LOGGER.info("Size plot saved as %s.png%s", output_prefix, " and .svg" if save_svg else "")

