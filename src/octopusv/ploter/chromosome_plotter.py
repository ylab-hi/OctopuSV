# chromosome_plotter.py

import logging

import matplotlib.pyplot as plt
import numpy as np

logging.basicConfig(level=logging.DEBUG)


class ChromosomePlotter:
    def __init__(self, input_file):
        self.input_file = input_file
        # GRCh38 chromosome lengths in base pairs
        self.chromosome_lengths = {
            "chr1": 248956422,
            "chr2": 242193529,
            "chr3": 198295559,
            "chr4": 190214555,
            "chr5": 181538259,
            "chr6": 170805979,
            "chr7": 159345973,
            "chr8": 145138636,
            "chr9": 138394717,
            "chr10": 133797422,
            "chr11": 135086622,
            "chr12": 133275309,
            "chr13": 114364328,
            "chr14": 107043718,
            "chr15": 101991189,
            "chr16": 90338345,
            "chr17": 83257441,
            "chr18": 80373285,
            "chr19": 58617616,
            "chr20": 64444167,
            "chr21": 46709983,
            "chr22": 50818468,
            "chrX": 156040895,
            "chrY": 57227415,
        }
        self.data = self.parse_data()

    def parse_data(self):
        """Parse the statistics file and extract chromosome data."""
        chromosome_order = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]
        chromosome_data = {chrom: {"count": 0, "density": 0} for chrom in chromosome_order}

        parsing_chromosomes = False
        with open(self.input_file) as f:
            for line in f:
                if "Chromosome Distribution" in line:
                    parsing_chromosomes = True
                    continue
                if parsing_chromosomes:
                    if line.strip() == "":
                        break
                    parts = line.strip().split("=")
                    if len(parts) == 2:
                        chrom = parts[0].strip()
                        if chrom in chromosome_data:
                            try:
                                count = int(parts[1].split()[0])
                                # Calculate density per Mb
                                length_mb = self.chromosome_lengths[chrom] / 1_000_000
                                density = count / length_mb
                                chromosome_data[chrom] = {
                                    "count": count,
                                    "density": round(density, 1),  # Round to 1 decimal place
                                }
                            except (ValueError, ZeroDivisionError) as e:
                                logging.warning(f"Error processing chromosome {chrom}: {e}")
                                chromosome_data[chrom] = {"count": 0, "density": 0}

        logging.debug(f"Parsed chromosome data: {chromosome_data}")
        return chromosome_data

    def plot(self, output_prefix, *, save_svg=True):
        """Create and save the chromosome distribution plot."""
        if not self.data:
            logging.error("No data to plot")
            return

        # Prepare data
        chromosomes = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]
        x = np.arange(len(chromosomes))  # x coordinates for the bars
        raw_counts = [self.data[chrom]["count"] for chrom in chromosomes]
        densities = [self.data[chrom]["density"] for chrom in chromosomes]

        # Create figure with two subplots
        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(20, 16), height_ratios=[1, 1])
        fig.suptitle("Structural Variant Distribution Across Chromosomes", fontsize=20, y=0.95, fontweight="bold")

        # Set the style for both plots
        for ax in [ax1, ax2]:
            ax.grid(True, axis="y", linestyle="--", alpha=0.3)
            ax.set_axisbelow(True)
            ax.spines["top"].set_visible(False)
            ax.spines["right"].set_visible(False)
            for spine in ax.spines.values():
                spine.set_color("#cccccc")
                spine.set_linewidth(0.8)
            ax.set_xticks(x)
            ax.tick_params(axis="both", which="major", labelsize=12)

        # Plot 1: Raw counts with refined color
        bars1 = ax1.bar(x, raw_counts, color="#4a90e2", width=0.8, edgecolor="white", linewidth=1)
        ax1.set_title("Raw SV Counts", fontsize=16, pad=20, fontweight="bold")
        ax1.set_ylabel("Number of SVs", fontsize=14, fontweight="bold")
        ax1.set_xticklabels(chromosomes, rotation=0)

        # Plot 2: Normalized densities with complementary color
        bars2 = ax2.bar(x, densities, color="#2ecc71", width=0.8, edgecolor="white", linewidth=1)
        ax2.set_title("Normalized SV Density", fontsize=16, pad=20, fontweight="bold")
        ax2.set_xlabel("Chromosome", fontsize=14, fontweight="bold")
        ax2.set_ylabel("SVs per Mb", fontsize=14, fontweight="bold")
        ax2.set_xticklabels(chromosomes, rotation=0)

        # Add value labels on bars
        def add_value_labels(axis, bars):
            for bar in bars:
                height = bar.get_height()
                axis.text(
                    bar.get_x() + bar.get_width() / 2.0,
                    height,
                    f"{height:,.1f}" if height % 1 else f"{int(height):,}",
                    ha="center",
                    va="bottom",
                    fontsize=10,
                    fontweight="bold",
                )

        add_value_labels(ax1, bars1)
        add_value_labels(ax2, bars2)

        # Adjust layout
        plt.tight_layout()

        # Save plots with high quality
        plt.savefig(f"{output_prefix}.png", dpi=300, bbox_inches="tight", facecolor="white")

        if save_svg:
            plt.savefig(f"{output_prefix}.svg", bbox_inches="tight", facecolor="white")

        plt.close()

        logging.info(f"Plot saved as {output_prefix}.png and {output_prefix}.svg")
