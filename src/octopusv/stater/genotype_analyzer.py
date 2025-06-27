from collections import Counter


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
                # Multiple callers case - each field is a separate caller
                gt = self._resolve_multi_caller_genotype(format_field, sample_data_fields, info_field)

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

    def _resolve_multi_caller_genotype(self, format_field, caller_segments, info_field):
        """Resolve genotype conflicts in caller mode using voting + AD tie-breaking + file order."""
        format_keys = format_field.split(":")

        if "GT" not in format_keys:
            return None

        gt_index = format_keys.index("GT")
        ad_index = format_keys.index("AD") if "AD" in format_keys else None

        # Collect genotypes and their supporting evidence
        genotype_data = []  # List of (genotype, ad_support, caller_index)

        for i, caller_segment in enumerate(caller_segments):
            sample_fields = caller_segment.split(":")
            if gt_index < len(sample_fields):
                gt = sample_fields[gt_index]

                # Get AD support if available
                ad_support = 0
                if ad_index is not None and ad_index < len(sample_fields):
                    ad_value = sample_fields[ad_index]
                    ad_support = self._extract_variant_support(ad_value)

                genotype_data.append((gt, ad_support, i))

        if not genotype_data:
            return None

        # Step 1: Voting - count occurrences of each genotype
        genotype_votes = [item[0] for item in genotype_data]
        vote_counts = Counter(genotype_votes)
        max_votes = max(vote_counts.values())
        tied_genotypes = [gt for gt, count in vote_counts.items() if count == max_votes]

        if len(tied_genotypes) == 1:
            # Clear winner from voting
            return tied_genotypes[0]

        # Step 2: Tie-breaking by AD support (variant supporting reads)
        # Prioritize genotypes with valid AD values (>= 0) over missing ones (-1)
        tied_genotype_ad = {}
        for gt, ad_support, caller_idx in genotype_data:
            if gt in tied_genotypes:
                if gt not in tied_genotype_ad or self._is_better_ad_support(ad_support, tied_genotype_ad[gt][0]):
                    tied_genotype_ad[gt] = (ad_support, caller_idx)

        # Find genotypes with the best AD support
        # First, separate genotypes with valid AD (>= 0) from those with missing AD (-1)
        valid_ad_genotypes = {gt: (ad, idx) for gt, (ad, idx) in tied_genotype_ad.items() if ad >= 0}
        missing_ad_genotypes = {gt: (ad, idx) for gt, (ad, idx) in tied_genotype_ad.items() if ad < 0}

        if valid_ad_genotypes:
            # If we have genotypes with valid AD values, use them
            max_ad = max(ad for ad, _ in valid_ad_genotypes.values())
            ad_winners = [gt for gt, (ad, _) in valid_ad_genotypes.items() if ad == max_ad]
        else:
            # All genotypes have missing AD, fall back to using all tied genotypes
            ad_winners = list(missing_ad_genotypes.keys())

        if len(ad_winners) == 1:
            # Clear winner from AD support
            return ad_winners[0]

        # Step 3: Final tie-breaking by SOURCES order (input file order)
        sources_order = self._extract_sources_order(info_field)

        if sources_order:
            # Try to map each tied genotype to its source position
            genotype_source_priority = {}

            for gt in ad_winners:
                # Get the caller index for this genotype
                caller_idx = tied_genotype_ad[gt][1]

                # Map caller index to source priority
                if caller_idx < len(sources_order):
                    # Find the position of this source in the original input order
                    source_name = sources_order[caller_idx]
                    genotype_source_priority[gt] = caller_idx
                else:
                    genotype_source_priority[gt] = 999  # Fallback for unknown sources

            # Return genotype with earliest source order
            best_gt = min(ad_winners, key=lambda gt: genotype_source_priority.get(gt, 999))
            return best_gt

        # Step 4: Final fallback - use first occurrence order in caller segments
        earliest_caller_idx = min(tied_genotype_ad[gt][1] for gt in ad_winners)
        for gt in ad_winners:
            if tied_genotype_ad[gt][1] == earliest_caller_idx:
                return gt

        # Ultimate fallback
        return ad_winners[0]

    def _extract_sources_order(self, info_field):
        """Extract SOURCES order from INFO field."""
        try:
            for item in info_field.split(";"):
                if item.startswith("SOURCES="):
                    sources = item.split("=", 1)[1].split(",")
                    return [s.strip() for s in sources]
            return []
        except:
            return []

    def _is_better_ad_support(self, new_ad, current_ad):
        """Determine if new AD support is better than current AD support.
        Valid AD values (>= 0) are always better than missing AD (-1).
        Among valid values, higher is better.
        """
        # If current AD is missing (-1) and new AD is valid (>= 0), new is better
        if current_ad < 0 and new_ad >= 0:
            return True
        # If new AD is missing (-1) and current AD is valid (>= 0), current is better
        if new_ad < 0 and current_ad >= 0:
            return False
        # If both are valid or both are missing, higher value is better
        return new_ad > current_ad

    def _extract_variant_support(self, ad_value):
        """Extract variant supporting read count from AD field (second number).
        Returns -1 for missing/empty AD values to deprioritize them.
        """
        try:
            # Handle empty or missing AD values
            if not ad_value or ad_value == "." or ad_value == "":
                return -1  # Missing AD should be deprioritized

            if "," in ad_value:
                ad_parts = ad_value.split(",")
                if len(ad_parts) >= 2:
                    second_part = ad_parts[1].strip()
                    # Check if second part is empty or "."
                    if second_part == "." or second_part == "":
                        return -1  # Missing variant support should be deprioritized
                    return int(second_part)  # Second number (variant support)
                else:
                    # Only one part, check if it's empty
                    first_part = ad_parts[0].strip()
                    if first_part == "." or first_part == "":
                        return -1
                    return int(first_part)
            elif ad_value != "." and ad_value != "":
                # Handle single value case
                return int(ad_value)
            else:
                return -1  # Empty or "." should be deprioritized
        except (ValueError, IndexError):
            return -1  # Invalid AD values should be deprioritized