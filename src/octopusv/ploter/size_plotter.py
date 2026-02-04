import logging

import matplotlib.pyplot as plt
import seaborn as sns

logging.basicConfig(level=logging.DEBUG)


class SizePlotter:
    def __init__(self, input_file):
        self.input_file = input_file
        self.data = self.parse_data()

    def parse_data(self):
        """Parse the statistics file and extract SV size distribution."""
        size_distribution = {}
        parsing_distribution = False
        with open(self.input_file) as f:
            for line in f:
                if "Size distribution" in line:
                    parsing_distribution = True
                    continue
                if parsing_distribution:
                    if line.strip() == "":
                        break
                    parts = line.strip().split("=")
                    if len(parts) == 2:
                        size_range = parts[0].strip()
                        count = int(parts[1].strip())
                        size_distribution[size_range] = count
        return size_distribution

    def plot(self, output_prefix, *, save_svg=True):
        """Create and save the SV size distribution plot."""
        if not self.data:
            logging.error("No data to plot")
            return

        # Create figure with specified size
        plt.figure(figsize=(12, 7))

        # Set style with refined grid
        sns.set_style("whitegrid", {"grid.linestyle": "--", "grid.alpha": 0.3})

        # Define consistent size order
        size_order = ["0-50 bp", "51-100 bp", "101-500 bp", "501-1 kb", "1 kb-10 kb", ">10 kb"]
        sizes = [self.data.get(size, 0) for size in size_order]

        # Create custom color gradient
        colors = sns.color_palette("RdYlBu_r", n_colors=len(size_order))

        # Create bar plot with refined styling
        bars = plt.bar(size_order, sizes, color=colors, width=0.7, edgecolor="white", linewidth=1.5)

        # Add value labels on bars
        for bar in bars:
            height = bar.get_height()
            plt.text(
                bar.get_x() + bar.get_width() / 2.0,
                height,
                f"{int(height):,}",
                ha="center",
                va="bottom",
                fontsize=10,
                fontweight="bold",
            )

        # Customize labels and title
        plt.xlabel("SV Size Range", fontsize=12, labelpad=10, fontweight="bold")
        plt.ylabel("Count", fontsize=12, labelpad=10, fontweight="bold")
        plt.title("Structural Variant Size Distribution", fontsize=14, pad=20, fontweight="bold")

        # Adjust tick labels
        plt.xticks(rotation=25, ha="right", fontsize=10)
        plt.yticks(fontsize=10)

        # Add subtle border
        for spine in plt.gca().spines.values():
            spine.set_visible(True)
            spine.set_color("#cccccc")
            spine.set_linewidth(0.8)

        # Adjust layout and save
        plt.tight_layout()
        plt.savefig(f"{output_prefix}.png", dpi=300, bbox_inches="tight")

        if save_svg:
            plt.savefig(f"{output_prefix}.svg", format="svg", bbox_inches="tight")

        plt.close()
