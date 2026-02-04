from collections import Counter


class TypeAnalyzer:
    def __init__(self, input_file):
        self.input_file = input_file

    def analyze(self):
        sv_types = Counter()
        with open(self.input_file) as f:
            for line in f:
                if line.startswith("#"):
                    continue
                fields = line.strip().split("\t")
                info = dict(item.split("=") for item in fields[7].split(";"))
                sv_type = info.get("SVTYPE", "Unknown")
                sv_types[sv_type] += 1

        total = sum(sv_types.values())
        result = "SV Type Analysis:\n"
        for sv_type, count in sv_types.items():
            percentage = (count / total) * 100
            result += f"{sv_type}: {count} ({percentage:.2f}%)\n"
        return result
