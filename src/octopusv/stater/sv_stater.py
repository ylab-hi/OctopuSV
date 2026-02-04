from pathlib import Path

from .chromosome_analyzer import ChromosomeAnalyzer
from .genotype_analyzer import GenotypeAnalyzer
from .qc_analyzer import QCAnalyzer
from .size_analyzer import SizeAnalyzer
from .type_analyzer import TypeAnalyzer


class SVStater:
    def __init__(self, input_file, min_size=50, max_size=None):
        self.input_file = input_file
        self.min_size = min_size
        self.max_size = max_size
        self.results = {}

    def analyze(self):
        self.results["type"] = TypeAnalyzer(self.input_file).analyze()
        self.results["size"] = SizeAnalyzer(self.input_file, self.min_size, self.max_size).analyze()
        self.results["chromosome"] = ChromosomeAnalyzer(self.input_file).analyze()
        self.results["qc"] = QCAnalyzer(self.input_file).analyze()
        self.results["genotype"] = GenotypeAnalyzer(self.input_file).analyze()

    def write_results(self, output_file):
        with Path.open(output_file, "w") as f:
            f.write("OctopusV report\n")
            f.write("-" * 40 + "\n\n")

            f.write(">>>>>>> Input\n\n")
            f.write(f"input file = {self.input_file}\n")
            f.write(f"output file = {output_file}\n\n")

            for category, result in self.results.items():
                f.write(f">>>>>>> {category.capitalize()} Analysis\n\n")
                f.write(self.format_result(result))
                f.write("\n")

    def format_result(self, result):
        formatted = ""
        for line in result.split("\n"):
            if ":" in line:
                key, value = line.split(":", 1)
                formatted += f"{key.strip():30} = {value.strip()}\n"
            else:
                formatted += line + "\n"
        return formatted

    def export_html(self, output_stat):
        """Export the analysis results to an HTML file."""
        result = {}
        result["input_file"] = self.input_file
        result["output_file"] = output_stat

        # Parse SV types
        result['sv_types'] = {}
        for line in self.results["type"].strip().split("\n")[1:]:
            if ":" in line:
                key, value = line.split(":", 1)
                parts = value.strip().split()
                if len(parts) >= 2 and "(" in parts[1]:
                    try:
                        count = int(parts[0])
                        percent = float(parts[1][1:-2])  # Remove percentage sign and parentheses
                        result['sv_types'][key.strip()] = (count, percent)
                    except (ValueError, IndexError):
                        continue

        # Parse size distribution
        lines = self.results["size"].strip().split("\n")
        try:
            result["total_svs"] = int(lines[0].split(":")[1].strip())
            result["min_size"] = int(lines[1].split(":")[1].strip().split()[0])
            result["max_size"] = int(lines[2].split(":")[1].strip().split()[0])
            result["mean_size"] = float(lines[3].split(":")[1].strip().split()[0])
            result["median_size"] = float(lines[4].split(":")[1].strip().split()[0])
            result['std_dev'] = float(lines[5].split(":")[1].strip().split()[0])
        except (IndexError, ValueError):
            # Set default values
            result["total_svs"] = 0
            result["min_size"] = self.min_size
            result["max_size"] = self.max_size or 0
            result["mean_size"] = 0.0
            result["median_size"] = 0.0
            result['std_dev'] = 0.0

        # Parse size distribution details
        result['size_distribution'] = {}
        for line in lines[8:]:
            if ":" in line:
                try:
                    key, value = line.split(":", 1)
                    result["size_distribution"][key.strip()] = int(value.strip())
                except (ValueError, IndexError):
                    continue

        # Parse QC results
        lines = self.results['qc'].strip().split("\n")

        # Safely get average quality value
        avg_qual_idx = -1
        for i, line in enumerate(lines):
            if "Average QUAL:" in line:
                avg_qual_idx = i
                break

        if avg_qual_idx >= 0:
            try:
                result['avg_qual'] = float(lines[avg_qual_idx].split(":")[1].strip())
            except (ValueError, IndexError):
                result['avg_qual'] = 0.0
        else:
            result['avg_qual'] = 0.0

        # Initialize filter status
        result["filter_status"] = {
            "PASS": ("0", 0.0),
            "hom_ref": ("0", 0.0),
            "not_fully_covered": ("0", 0.0)
        }

        # Find filter status section
        filter_idx = -1
        for i, line in enumerate(lines):
            if "Filter Status:" in line:
                filter_idx = i
                break

        # Parse filter status information
        if filter_idx >= 0:
            filter_lines = []
            for i in range(filter_idx + 1, len(lines)):
                if lines[i].strip() and ":" in lines[i]:
                    filter_lines.append(lines[i])
                elif lines[i].strip() == "":  # Empty line indicates end of section
                    break

            # Process filter status
            for line in filter_lines:
                try:
                    status, value = line.split(":", 1)
                    status = status.strip()
                    value_parts = value.strip().split()

                    if len(value_parts) >= 1:
                        count = value_parts[0]
                        percent = 0.0

                        if len(value_parts) >= 2 and "(" in value_parts[1]:
                            try:
                                percent = float(value_parts[1][1:-2])  # Remove parentheses and percentage sign
                            except ValueError:
                                percent = 0.0

                        # Update or add status
                        result["filter_status"][status] = (count, percent)
                except (ValueError, IndexError):
                    continue

        # Safely get average read support
        read_support_idx = -1
        for i, line in enumerate(lines):
            if "Average Read Support:" in line:
                read_support_idx = i
                break

        if read_support_idx >= 0:
            try:
                value = lines[read_support_idx].split(":")[1].strip()
                # Check if it contains percentage
                if " " in value and "(" in value:
                    value = value.split()[0]  # Only take the first part (number)
                result['avg_read_support'] = float(value)
            except (ValueError, IndexError):
                result['avg_read_support'] = 0.0
        else:
            result['avg_read_support'] = 0.0

        # Parse genotype distribution (enhanced for multi-sample support)
        lines = self.results['genotype'].strip().split("\n")
        result["genotype_dist"] = {}
        result["sample_genotypes"] = {}  # For per-sample genotypes
        result["population_genotypes"] = {}  # For population-level genotypes

        current_sample = None
        is_population = False
        is_in_genotype_section = False

        for line in lines:
            line = line.strip()

            # Skip empty lines
            if not line:
                continue

            # Check for main header
            if line == "Genotype Distribution:":
                is_in_genotype_section = True
                continue

            if not is_in_genotype_section:
                continue

            # Check for sample-specific genotype sections
            if " Genotypes:" in line and not line.startswith("Overall"):
                current_sample = line.split(" Genotypes:")[0].strip()
                result["sample_genotypes"][current_sample] = {}
                is_population = False
                continue
            elif line.startswith("Overall:"):
                current_sample = None
                is_population = True
                continue

            # Parse genotype data lines
            if ":" in line:
                try:
                    key, value = line.split(":", 1)
                    key = key.strip()
                    value = value.strip()

                    # Remove leading spaces for indented lines
                    if key.startswith("  "):
                        key = key[2:]

                    # Parse value (count and percentage)
                    if " " in value and "(" in value:
                        parts = value.split()
                        count = int(parts[0])
                        percent_str = parts[1]
                        if percent_str.startswith("(") and percent_str.endswith("%)"):
                            percent = float(percent_str[1:-2])  # Remove (%) wrapper
                        else:
                            percent = 0.0
                        genotype_data = (count, percent)
                    else:
                        try:
                            genotype_data = (int(value), 0.0)
                        except ValueError:
                            continue

                    # Store in appropriate location
                    if current_sample:
                        # Sample-specific genotype
                        result["sample_genotypes"][current_sample][key] = genotype_data
                    elif is_population:
                        # Population-level genotype
                        result["population_genotypes"][key] = genotype_data
                    else:
                        # Legacy single-sample format (caller mode)
                        result["genotype_dist"][key] = genotype_data

                except (ValueError, IndexError):
                    continue

        return result