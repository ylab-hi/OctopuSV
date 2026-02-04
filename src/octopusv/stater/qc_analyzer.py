import statistics
from collections import defaultdict


class QCAnalyzer:
    def __init__(self, input_file):
        self.input_file = input_file
        self.default_qual = 0  # Default quality score when missing

    def _safe_float_convert(self, value, default=0):
        """Safely convert a value to float, returning default if conversion fails."""
        try:
            if value in (".", ""):
                return default
            return float(value)
        except (ValueError, TypeError):
            return default

    def _safe_int_convert(self, value, default=0):
        """Safely convert a value to integer, returning default if conversion fails."""
        try:
            if value in (".", ""):
                return default
            return int(value)
        except (ValueError, TypeError):
            return default

    def _calculate_statistics(self, values):
        """Calculate statistics safely handling empty lists."""
        if not values:
            return {"mean": 0, "median": 0, "min": 0, "max": 0, "q1": 0, "q3": 0}

        try:
            return {
                "mean": statistics.mean(values),
                "median": statistics.median(values),
                "min": min(values),
                "max": max(values),
                "q1": statistics.quantiles(values, n=4)[0] if len(values) >= 4 else values[0],
                "q3": statistics.quantiles(values, n=4)[2] if len(values) >= 4 else values[-1],
            }
        except statistics.StatisticsError:
            return {"mean": 0, "median": 0, "min": 0, "max": 0, "q1": 0, "q3": 0}

    def _parse_info_field(self, info_str):
        """Safely parse the INFO field."""
        info_dict = {}
        try:
            for item in info_str.split(";"):
                if "=" in item:
                    key, value = item.split("=", 1)
                    info_dict[key] = value
                else:
                    info_dict[item] = True
        except Exception:
            pass
        return info_dict

    def analyze(self):
        qual_scores = []
        filter_status = defaultdict(int)
        support_counts = []

        try:
            with open(self.input_file) as f:
                for line in f:
                    if line.startswith("#"):
                        continue

                    try:
                        fields = line.strip().split("\t")
                        if len(fields) < 8:  # Ensure minimum required fields
                            continue

                        # Handle quality score
                        qual = self._safe_float_convert(fields[5])
                        if qual != 0:  # Only append non-zero values
                            qual_scores.append(qual)

                        # Handle filter status
                        filter_val = fields[6] if len(fields) > 6 else "UNKNOWN"
                        filter_status[filter_val] += 1

                        # Handle support count from INFO field
                        if len(fields) > 7:
                            info = self._parse_info_field(fields[7])
                            support = self._safe_int_convert(info.get("SUPPORT", "0"))
                            if support > 0:  # Only append non-zero values
                                support_counts.append(support)

                    except IndexError:
                        continue  # Skip malformed lines

        except FileNotFoundError:
            return "Error: Input file not found"
        except Exception as e:
            return f"Error: Failed to analyze file: {e!s}"

        # Calculate statistics
        qual_stats = self._calculate_statistics(qual_scores)
        support_stats = self._calculate_statistics(support_counts)

        # Format results
        result = "Quality Score (QUAL) Statistics:\n"
        result += f"Number of variants with QUAL: {len(qual_scores)}\n"
        result += f"Average QUAL: {qual_stats['mean']:.2f}\n"
        result += f"Median QUAL: {qual_stats['median']:.2f}\n"
        result += f"Min QUAL: {qual_stats['min']:.2f}\n"
        result += f"Max QUAL: {qual_stats['max']:.2f}\n"
        result += f"Q1 QUAL: {qual_stats['q1']:.2f}\n"
        result += f"Q3 QUAL: {qual_stats['q3']:.2f}\n\n"

        result += "Filter Status:\n"
        total_variants = sum(filter_status.values())
        for status, count in sorted(filter_status.items()):
            percentage = (count / total_variants * 100) if total_variants > 0 else 0
            result += f"{status}: {count} ({percentage:.2f}%)\n"

        result += "\nRead Support Statistics:\n"
        result += f"Number of variants with support info: {len(support_counts)}\n"
        result += f"Average Read Support: {support_stats['mean']:.2f}\n"
        result += f"Median Read Support: {support_stats['median']:.2f}\n"
        result += f"Min Read Support: {support_stats['min']:.2f}\n"
        result += f"Max Read Support: {support_stats['max']:.2f}\n"
        result += f"Q1 Read Support: {support_stats['q1']:.2f}\n"
        result += f"Q3 Read Support: {support_stats['q3']:.2f}\n"

        return result
