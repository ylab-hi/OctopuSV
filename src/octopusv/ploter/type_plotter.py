import collections
import logging

import matplotlib.pyplot as plt

logging.basicConfig(level=logging.DEBUG)

# Define consistent color scheme for SV types
COLOR_SCHEME = {
    "TRA": "#FF6B6B",  # Coral red
    "INV": "#4ECDC4",  # Turquoise
    "DUP": "#45B7D1",  # Sky blue
    "INS": "#96CEB4",  # Mint green
    "DEL": "#6C5B7B",  # Deep purple
}


class TypePlotter:
    def __init__(self, input_file):
        self.input_file = input_file
        self.sv_order = ["TRA", "INV", "DUP", "INS", "DEL"]
        self.color_map = COLOR_SCHEME
        self.data = self.parse_data()

    def parse_data(self):
        """Parse the statistics file and extract SV type information."""
        sv_types = {}
        parsing_types = False
        with open(self.input_file) as f:
            for line in f:
                if "SV Type Analysis" in line:
                    parsing_types = True
                    continue
                if parsing_types:
                    if line.strip() == "":
                        break
                    parts = line.strip().split("=")
                    if len(parts) == 2:
                        sv_type = parts[0].strip()
                        count_percentage = parts[1].strip().split()
                        if len(count_percentage) >= 2:
                            count = int(count_percentage[0])
                            percentage = float(count_percentage[1].strip("()%"))
                            sv_types[sv_type] = (count, percentage)
        return sv_types

    def _sort_data(self):
        """Sort data according to predefined order."""
        ordered_data = collections.OrderedDict()
        for sv_type in self.sv_order:
            if sv_type in self.data:
                ordered_data[sv_type] = self.data[sv_type]
        for sv_type, values in self.data.items():
            if sv_type not in self.sv_order:
                ordered_data[sv_type] = values
        return ordered_data

    def plot(self, output_prefix, *, save_svg=True):
        """Create and save the SV type distribution plot."""
        if not self.data:
            logging.error("No data to plot")
            return

        ordered_data = self._sort_data()

        # Set up the figure
        plt.figure(figsize=(10, 8))
        plt.rcParams["font.family"] = "sans-serif"
        plt.rcParams["font.sans-serif"] = ["Arial"]
        plt.rcParams["svg.fonttype"] = "none"  # Ensure editable text in SVG

        # Prepare data
        types = list(ordered_data.keys())
        sizes = [count for count, _ in ordered_data.values()]
        colors = [self.color_map.get(sv_type, "#7f8c8d") for sv_type in types]

        # Create pie chart with refined styling
        explode = [0.01] * len(types)  # Slight separation between segments
        wedges, texts, autotexts = plt.pie(
            sizes,
            explode=explode,
            labels=types,
            colors=colors,
            autopct="%1.1f%%",
            pctdistance=0.85,
            startangle=90,
            wedgeprops={"width": 0.5, "edgecolor": "white", "linewidth": 2},
        )

        # Style percentage labels
        plt.setp(autotexts, size=9, weight="bold", color="white")
        plt.setp(texts, size=10, weight="bold")

        # Add center circle for donut effect
        centre_circle = plt.Circle((0, 0), 0.70, fc="white", ec="#e0e0e0")
        fig = plt.gcf()
        fig.gca().add_artist(centre_circle)

        # Configure legend
        legend_labels = [f"{sv_type}: {ordered_data[sv_type][0]:,}" for sv_type in types]
        plt.legend(
            wedges,
            legend_labels,
            title="SV Types",
            loc="center left",
            bbox_to_anchor=(1, 0, 0.5, 1),
            frameon=False,
            title_fontsize=12,
            fontsize=10,
        )

        # Add title and ensure proportional scaling
        plt.title("Structural Variant Type Distribution", pad=20, fontsize=14, fontweight="bold")
        plt.axis("equal")

        # Save plots
        plt.savefig(f"{output_prefix}.png", dpi=300, bbox_inches="tight", facecolor="white")

        if save_svg:
            plt.savefig(f"{output_prefix}.svg", bbox_inches="tight", facecolor="white")
        plt.close()
