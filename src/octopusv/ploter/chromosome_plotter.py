"""Chromosome distribution plotter.

The new stat pipeline passes a structured chromosome stats dict directly:

    {
        "counts": {"1": 3809, "2": 4038, ...},
        "density_per_mb": {"1": 15.30, "2": 16.60, ...},
        "unmatched_contigs": [...]
    }

This plotter no longer hardcodes genome lengths and no longer parses stat.txt.
For backward compatibility, it can still accept a path to a JSON stat output or
a legacy text stat file, but the preferred path is dict input from SVStater.
"""

from __future__ import annotations

import json
import logging
import re
from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt
import numpy as np

LOGGER = logging.getLogger(__name__)


def _natural_chrom_key(chrom: str) -> tuple[int, int | str]:
    """Natural chromosome order: 1..22, X, Y, MT, then other contigs."""
    raw = str(chrom)
    name = raw[3:] if raw.lower().startswith("chr") else raw

    if name.isdigit():
        value = int(name)
        if 1 <= value <= 22:
            return (0, value)

    if name == "X":
        return (0, 23)
    if name == "Y":
        return (0, 24)
    if name in {"M", "MT"}:
        return (0, 25)

    return (1, raw)


def _read_json(path: Path) -> dict[str, Any]:
    with path.open() as handle:
        return json.load(handle)


class ChromosomePlotter:
    """Plot chromosome-level SV counts and optional density."""

    def __init__(self, stats_or_file: dict[str, Any] | str | Path):
        self.stats_or_file = stats_or_file
        self.data = self._load_data(stats_or_file)

    def _load_data(self, stats_or_file: dict[str, Any] | str | Path) -> dict[str, Any]:
        if isinstance(stats_or_file, dict):
            return self._from_dict(stats_or_file)

        path = Path(stats_or_file)
        if path.suffix.lower() == ".json":
            return self._from_dict(_read_json(path))

        return self._from_legacy_text(path)

    def _from_dict(self, data: dict[str, Any]) -> dict[str, Any]:
        """Normalize either full stat JSON or chromosome-only stats."""
        chrom_stats = data.get("chromosome", data)

        counts = chrom_stats.get("counts", {}) or {}
        density = chrom_stats.get("density_per_mb", {}) or {}
        unmatched = chrom_stats.get("unmatched_contigs", []) or []
        length_source = chrom_stats.get("length_source")

        return {
            "counts": {str(k): int(v) for k, v in counts.items()},
            "density_per_mb": {
                str(k): float(v) for k, v in density.items() if v is not None
            },
            "unmatched_contigs": [str(x) for x in unmatched],
            "length_source": length_source,
        }

    def _from_legacy_text(self, path: Path) -> dict[str, Any]:
        """Best-effort parser for old stat.txt files.

        This is only for backward compatibility. The preferred input is a dict.
        """
        counts: dict[str, int] = {}
        density: dict[str, float] = {}

        in_section = False
        pattern = re.compile(
            r"^\s*(?P<chrom>[^=]+?)\s*=\s*"
            r"(?P<count>\d+)\s+SVs"
            r"(?:\s+\((?P<density>[0-9.]+)\s+SVs/Mb.*\))?"
        )

        with path.open() as handle:
            for line in handle:
                if "Chromosome Distribution" in line:
                    in_section = True
                    continue
                if in_section and not line.strip():
                    break
                if not in_section:
                    continue

                match = pattern.match(line)
                if not match:
                    continue

                chrom = match.group("chrom").strip()
                counts[chrom] = int(match.group("count"))

                if match.group("density") is not None:
                    density[chrom] = float(match.group("density"))

        return {
            "counts": counts,
            "density_per_mb": density,
            "unmatched_contigs": [],
            "length_source": None,
        }

    def _ordered_counts(self) -> list[tuple[str, int]]:
        return sorted(self.data["counts"].items(), key=lambda kv: _natural_chrom_key(kv[0]))

    def _ordered_density(self) -> list[tuple[str, float]]:
        density = self.data.get("density_per_mb") or {}
        return sorted(density.items(), key=lambda kv: _natural_chrom_key(kv[0]))

    def plot(self, output_prefix: str | Path, *, save_svg: bool = True) -> None:
        """Create and save chromosome count/density plot."""
        ordered_counts = self._ordered_counts()
        if not ordered_counts:
            LOGGER.error("No chromosome count data to plot.")
            return

        ordered_density = self._ordered_density()
        has_density = bool(ordered_density)

        if has_density:
            fig, axes = plt.subplots(2, 1, figsize=(22, 16), height_ratios=[1, 1])
            ax_count, ax_density = axes
        else:
            fig, ax_count = plt.subplots(1, 1, figsize=(22, 8))
            ax_density = None

        fig.suptitle(
            "Structural Variant Distribution Across Chromosomes",
            fontsize=20,
            y=0.98,
            fontweight="bold",
        )

        # Raw counts
        chroms = [chrom for chrom, _ in ordered_counts]
        counts = [count for _, count in ordered_counts]
        x = np.arange(len(chroms))

        bars = ax_count.bar(x, counts, width=0.8, edgecolor="white", linewidth=1)
        ax_count.set_title("Raw SV Counts", fontsize=16, pad=18, fontweight="bold")
        ax_count.set_ylabel("Number of SVs", fontsize=14, fontweight="bold")
        ax_count.set_xticks(x)
        ax_count.set_xticklabels(chroms, rotation=45, ha="right")
        ax_count.grid(True, axis="y", linestyle="--", alpha=0.3)
        ax_count.set_axisbelow(True)
        ax_count.spines["top"].set_visible(False)
        ax_count.spines["right"].set_visible(False)

        if len(chroms) <= 35:
            self._add_value_labels(ax_count, bars, integer=True)

        # Density, only where density was actually computed.
        if ax_density is not None:
            density_chroms = [chrom for chrom, _ in ordered_density]
            densities = [value for _, value in ordered_density]
            x2 = np.arange(len(density_chroms))

            bars2 = ax_density.bar(x2, densities, width=0.8, edgecolor="white", linewidth=1)
            ax_density.set_title(
                "Normalized SV Density",
                fontsize=16,
                pad=18,
                fontweight="bold",
            )
            ax_density.set_xlabel("Chromosome", fontsize=14, fontweight="bold")
            ax_density.set_ylabel("SVs per Mb", fontsize=14, fontweight="bold")
            ax_density.set_xticks(x2)
            ax_density.set_xticklabels(density_chroms, rotation=45, ha="right")
            ax_density.grid(True, axis="y", linestyle="--", alpha=0.3)
            ax_density.set_axisbelow(True)
            ax_density.spines["top"].set_visible(False)
            ax_density.spines["right"].set_visible(False)

            if len(density_chroms) <= 35:
                self._add_value_labels(ax_density, bars2, integer=False)

        plt.tight_layout()
        output_prefix = str(output_prefix)
        plt.savefig(f"{output_prefix}.png", dpi=300, bbox_inches="tight", facecolor="white")

        if save_svg:
            plt.savefig(f"{output_prefix}.svg", bbox_inches="tight", facecolor="white")

        plt.close(fig)
        LOGGER.info("Chromosome plot saved as %s.png%s", output_prefix, " and .svg" if save_svg else "")

    @staticmethod
    def _add_value_labels(axis, bars, *, integer: bool) -> None:
        for bar in bars:
            height = bar.get_height()
            if height == 0:
                continue

            label = f"{int(height):,}" if integer else f"{height:.2f}"
            axis.text(
                bar.get_x() + bar.get_width() / 2.0,
                height,
                label,
                ha="center",
                va="bottom",
                fontsize=9,
                fontweight="bold",
            )
