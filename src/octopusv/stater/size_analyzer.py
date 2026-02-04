import statistics
from collections import defaultdict


class SizeAnalyzer:
    def __init__(self, input_file, min_size=50, max_size=None):
        self.input_file = input_file
        self.min_size = min_size
        self.max_size = max_size

    def analyze(self):
        sizes = []
        size_distribution = defaultdict(int)

        with open(self.input_file) as f:
            for line in f:
                if line.startswith("#"):
                    continue
                fields = line.strip().split("\t")
                info = dict(item.split("=") for item in fields[7].split(";"))

                if "SVLEN" in info and info["SVLEN"] != ".":
                    size = abs(int(info["SVLEN"]))
                else:
                    try:
                        size = abs(int(info.get("END", fields[1])) - int(fields[1]))
                    except ValueError:
                        continue

                if size < self.min_size or (self.max_size and size > self.max_size):
                    continue

                sizes.append(size)
                if size <= 50:
                    size_distribution["0-50 bp"] += 1
                elif size <= 100:
                    size_distribution["51-100 bp"] += 1
                elif size <= 500:
                    size_distribution["101-500 bp"] += 1
                elif size <= 1000:
                    size_distribution["501-1 kb"] += 1
                elif size <= 10000:
                    size_distribution["1 kb-10 kb"] += 1
                else:
                    size_distribution[">10 kb"] += 1

        result = f"Total SVs analyzed: {len(sizes)}\n"
        result += f"Minimum size: {min(sizes)} bp\n"
        result += f"Maximum size: {max(sizes)} bp\n"
        result += f"Mean size: {statistics.mean(sizes):.2f} bp\n"
        result += f"Median size: {statistics.median(sizes)} bp\n"
        result += f"Standard deviation: {statistics.stdev(sizes):.2f} bp\n\n"
        result += "Size distribution:\n"
        for range, count in size_distribution.items():
            result += f"{range}: {count}\n"

        return result
