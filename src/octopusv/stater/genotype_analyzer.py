from collections import Counter

from octopusv.utils.genotype_resolver import resolve_multi_caller_genotype


class GenotypeAnalyzer:
    def __init__(self, input_file):
        self.input_file = input_file

    def analyze(self):
        """Analyze genotype distribution with smart resolution for different modes."""

        with open(self.input_file) as f:
            header_line = None
            sample_names = []

            # Find header line and extract sample names
            for line in f:
                if line.startswith("#CHROM"):
                    header_line = line.strip()
                    header_fields = header_line.split("\t")
                    sample_names = header_fields[9:]  # Get sample names
                    break
                if line.startswith("#"):
                    continue

            # Reset file pointer to read data lines
            f.seek(0)

            # Determine mode based on number of sample columns
            if len(sample_names) == 1:
                # Caller mode: single SAMPLE column (may contain multiple caller data)
                return self._analyze_caller_mode(f)
            else:
                # Sample mode: multiple sample columns
                return self._analyze_sample_mode(f, sample_names)

    def _analyze_caller_mode(self, file_handle):
        """Analyze genotype distribution for caller mode."""
        genotypes = Counter()

        for line in file_handle:
            if line.startswith("#"):
                continue

            fields = line.strip().split("\t")
            if len(fields) < 10:
                continue

            format_field = fields[8]
            info_field = fields[7]

            # In caller mode, multiple caller data may be split across multiple fields
            # starting from field 9. We need to collect all sample data fields.
            sample_data_fields = fields[9:]

            if len(sample_data_fields) == 1:
                # Single caller case
                gt = self._resolve_single_caller_genotype(format_field, sample_data_fields[0])
            else:
                # Multiple callers case - resolve with the shared voting rule
                gt = resolve_multi_caller_genotype(format_field, sample_data_fields, info_field)

            if gt:
                genotypes[gt] += 1

        # Format results
        total = sum(genotypes.values())
        if total == 0:
            return "Genotype Distribution:\nNo valid genotypes found\n"

        result = "Genotype Distribution:\n"
        for gt, count in genotypes.items():
            percentage = (count / total) * 100
            result += f"{gt}: {count} ({percentage:.2f}%)\n"

        return result

    def _analyze_sample_mode(self, file_handle, sample_names):
        """Analyze genotype distribution for sample mode with per-sample and population stats."""
        # Initialize counters for each sample
        sample_genotypes = {name: Counter() for name in sample_names}

        for line in file_handle:
            if line.startswith("#"):
                continue

            fields = line.strip().split("\t")
            if len(fields) < 9 + len(sample_names):
                continue

            format_field = fields[8]
            sample_columns = fields[9:9 + len(sample_names)]

            # Extract genotype for each sample
            for i, sample_name in enumerate(sample_names):
                gt = self._resolve_single_caller_genotype(format_field, sample_columns[i])
                if gt:
                    sample_genotypes[sample_name][gt] += 1

        # Format results
        result = "Genotype Distribution:\n\n"

        # Per-sample statistics
        population_totals = Counter()
        for sample_name, genotypes in sample_genotypes.items():
            total = sum(genotypes.values())
            if total > 0:
                result += f"{sample_name} Genotypes:\n"
                for gt, count in genotypes.items():
                    percentage = (count / total) * 100
                    result += f"  {gt}: {count} ({percentage:.2f}%)\n"
                    population_totals[gt] += count
                result += "\n"

        # Population-level statistics
        pop_total = sum(population_totals.values())
        if pop_total > 0:
            result += "Overall:\n"
            for gt, count in population_totals.items():
                percentage = (count / pop_total) * 100
                result += f"  {gt}: {count} ({percentage:.2f}%)\n"

        return result

    def _resolve_single_caller_genotype(self, format_field, sample_data):
        """Resolve genotype for single caller data."""
        format_keys = format_field.split(":")

        if "GT" not in format_keys:
            return None

        gt_index = format_keys.index("GT")
        sample_fields = sample_data.split(":")

        if gt_index < len(sample_fields):
            return sample_fields[gt_index]

        return None
