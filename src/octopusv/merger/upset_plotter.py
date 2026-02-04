import os

import matplotlib.pyplot as plt
import numpy as np


class UpSetPlotter:
    def __init__(self, merged_events, input_files):
        """Initialize UpSetPlotter with merged events and input files."""
        self.input_files = [os.path.basename(f) for f in input_files]
        self.merged_events = merged_events
        self.intersections = self._get_intersections()
        self.file_totals = self._get_file_totals()
        self.short_input_files = self._shorten_file_names(self.input_files)

        # Color settings
        self.matrix_dot_color = "black"
        self.matrix_bg_dot_color = "#D3D3D3"
        self.matrix_line_color = "#404040"
        self.bar_color = "#4e79a7"
        self.intersection_bar_color = "#404040"
        self.alternate_row_color = "#f5f5f5"  # 非常浅的灰色作为交替行的背景色

    def _shorten_file_names(self, file_names, max_length=30):
        """Shorten file names if they exceed the maximum length."""
        shortened_names = []
        for name in file_names:
            shortened = name[:15] + "..." + name[-10:] if len(name) > max_length else name
            shortened_names.append(shortened)
        return shortened_names

    def _get_intersections(self):
        """Calculate all possible intersections and their sizes."""
        intersections = {}
        for event in self.merged_events:
            event_sources = {os.path.basename(s) for s in event.source_file.split(",")}
            source_key = tuple(sorted(event_sources))
            intersections[source_key] = intersections.get(source_key, 0) + 1
        return intersections

    def _get_file_totals(self):
        """Calculate the total number of events in each input file."""
        file_totals = {file: 0 for file in self.input_files}
        for event in self.merged_events:
            event_sources = {os.path.basename(s) for s in event.source_file.split(",")}
            for file in event_sources:
                if file in file_totals:
                    file_totals[file] += 1
        return file_totals

    def plot(self, output_file):
        """Create and save the UpSet plot."""
        # Sort input files by total size (from most to least)
        self.input_files.sort(key=lambda x: self.file_totals[x], reverse=True)
        self.short_input_files = [self._shorten_file_names([name])[0] for name in self.input_files]

        num_files = len(self.input_files)

        # Get intersection combinations sorted by size (from high to low)
        sorted_intersections = sorted(self.intersections.items(), key=lambda x: (x[1], len(x[0])), reverse=True)

        num_intersections = len(sorted_intersections)

        # Setup figure dimensions and spacing
        fig_width = max(10, num_intersections * 0.5)
        fig_height = max(8, num_files * 0.8)

        # Create figure with tighter spacing
        fig = plt.figure(figsize=(fig_width, fig_height))
        gs = fig.add_gridspec(
            nrows=3, ncols=2, height_ratios=[1, 0.05, 2], width_ratios=[1, 2.5], hspace=0.1, wspace=0.2
        )

        # Set up axes
        set_size_ax = fig.add_subplot(gs[2, 0])
        inter_size_ax = fig.add_subplot(gs[0, 1])
        matrix_ax = fig.add_subplot(gs[2, 1])

        # Plot set sizes (horizontal bars on left)
        set_sizes = [self.file_totals[file] for file in self.input_files]
        y_pos = np.arange(num_files)
        bars = set_size_ax.barh(y_pos, set_sizes, color=self.bar_color, height=0.6)

        # Remove y-axis ticks and labels
        set_size_ax.set_yticks(y_pos)
        set_size_ax.set_yticklabels([])
        set_size_ax.invert_yaxis()
        set_size_ax.set_xlabel("Set Size", fontsize=10, labelpad=10)
        set_size_ax.spines["top"].set_visible(False)
        set_size_ax.spines["right"].set_visible(False)

        # Add sample names to the right of the bars
        for bar, name in zip(bars, self.short_input_files, strict=False):
            width = bar.get_width()
            y = bar.get_y() + bar.get_height() / 2
            set_size_ax.text(width + max(set_sizes) * 0.02, y, name, va="center", ha="left", fontsize=10)

        # Adjust x-axis limits to make room for labels
        set_size_ax.set_xlim(0, max(set_sizes) * 1.2)

        # Plot intersection sizes (vertical bars on top)
        intersection_sizes = [count for _, count in sorted_intersections]
        x_pos = np.arange(len(sorted_intersections))
        bars = inter_size_ax.bar(x_pos, intersection_sizes, color=self.intersection_bar_color, width=0.6)

        # Add value labels to intersection bars
        for bar in bars:
            height = bar.get_height()
            inter_size_ax.text(
                bar.get_x() + bar.get_width() / 2,
                height + (max(intersection_sizes) * 0.02),
                f"{int(height)}",
                ha="center",
                va="bottom",
                fontsize=8,
                color="gray",
            )

        # Customize intersection size axis
        inter_size_ax.set_xticks([])
        inter_size_ax.set_ylabel("Intersection Size", fontsize=10, labelpad=10)
        inter_size_ax.spines["top"].set_visible(False)
        inter_size_ax.spines["right"].set_visible(False)

        # Create and plot matrix
        matrix = np.zeros((num_files, num_intersections))
        for i, (sources, _) in enumerate(sorted_intersections):
            for j, file in enumerate(self.input_files):
                if file in sources:
                    matrix[j, i] = 1

        # Add alternating row backgrounds
        for i in range(num_files):
            if i % 2 == 0:  # 偶数行添加背景
                matrix_ax.axhspan(i - 0.5, i + 0.5, color=self.alternate_row_color, zorder=0)

        # Plot background dots for all positions
        for i in range(num_files):
            for j in range(num_intersections):
                matrix_ax.scatter(j, i, color=self.matrix_bg_dot_color, s=60, zorder=1)

        # Plot black dots and connecting lines
        for j in range(num_intersections):
            indices = np.where(matrix[:, j])[0]
            # Plot black dots for connected positions
            for idx in indices:
                matrix_ax.scatter(j, idx, color=self.matrix_dot_color, s=60, zorder=3)
            # Draw connecting lines between dots
            if len(indices) > 1:
                matrix_ax.plot(
                    [j] * len(indices),
                    indices,
                    color=self.matrix_line_color,
                    linewidth=1,
                    solid_capstyle="round",
                    zorder=2,
                )

        # Customize matrix appearance
        matrix_ax.set_xlim(-0.5, num_intersections - 0.5)
        matrix_ax.set_ylim(num_files - 0.5, -0.5)
        matrix_ax.set_xticks([])
        matrix_ax.set_yticks([])
        matrix_ax.set_facecolor("white")

        # Remove matrix border
        for spine in matrix_ax.spines.values():
            spine.set_visible(False)

        # Adjust layout to make room for labels
        plt.subplots_adjust(left=0.1, right=0.9, wspace=0.3)

        # Save plot with white background
        plt.savefig(output_file, dpi=300, bbox_inches="tight", facecolor="white")
        plt.close()
