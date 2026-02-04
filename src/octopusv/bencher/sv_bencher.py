import logging
from pathlib import Path

from octopusv.utils.svcf_parser import SVCFFileEventCreator

from .bench_utils import calculate_metrics, write_summary, write_vcf


class SVBencher:
    """Benchmark structural variants using GIAB standards."""

    def __init__(
        self,
        truth_file: Path,
        call_file: Path,
        output_dir: Path,
        reference_distance: int = 500,
        sequence_similarity: float = 0.7,
        size_similarity: float = 0.7,
        reciprocal_overlap: float = 0.0,
        type_ignore: bool = False,
        size_min: int = 50,
        size_max: int = 50000,
        pass_only: bool = False,
        enable_sequence_comparison: bool = False,
    ):
        """Initialize the SV benchmarker with GIAB standard parameters."""
        self.truth_file = truth_file
        self.call_file = call_file
        self.output_dir = output_dir
        self.reference_distance = reference_distance
        self.sequence_similarity = sequence_similarity
        self.size_similarity = size_similarity
        self.reciprocal_overlap = reciprocal_overlap
        self.type_ignore = type_ignore
        self.size_min = size_min
        self.size_max = size_max
        self.pass_only = pass_only
        self.enable_sequence_comparison = enable_sequence_comparison

        self.truth_events = None
        self.call_events = None
        self.results = None

        self.logger = logging.getLogger(__name__)

    def run_benchmark(self):
        """Run the complete benchmarking process."""
        try:
            self._parse_files()
            self._compare_events()
            self._write_results()
        except Exception as e:
            self.logger.error(f"Benchmarking failed: {e!s}")
            raise

    def _parse_files(self):
        """Parse input VCF files."""
        self.logger.info("Parsing truth file...")
        truth_parser = SVCFFileEventCreator([str(self.truth_file)])
        truth_parser.parse()

        self.logger.info("Parsing call file...")
        call_parser = SVCFFileEventCreator([str(self.call_file)])
        call_parser.parse()

        self.truth_events = truth_parser.events
        self.call_events = call_parser.events

    def _filter_events(self, events: list) -> list:
        """Filter events based on size and FILTER criteria."""
        filtered = []
        for event in events:
            # Skip if not PASS and pass_only is True
            if self.pass_only and event.filter != "PASS":
                continue

            # Get event size
            try:
                if event.sv_type == "TRA":
                    size = 0  # TRA events don't have a meaningful size
                else:
                    size = abs(event.end_pos - event.start_pos)

                # Skip if outside size range (except for TRA)
                if event.sv_type != "TRA" and (size < self.size_min or size > self.size_max):
                    continue

                filtered.append(event)
            except AttributeError as e:
                self.logger.warning(f"Skipping malformed event: {e!s}")
                continue

        return filtered

    def _meets_matching_criteria(self, truth_event, call_event) -> bool:
        """Check if two events meet the matching criteria based on GIAB standards."""
        try:
            # Check SV type unless ignored
            if not self.type_ignore and truth_event.sv_type != call_event.sv_type:
                return False

            # Special handling for translocations
            if truth_event.sv_type == "TRA":
                return self._compare_tra_events(truth_event, call_event)

            # Check reference distance
            start_dist = abs(truth_event.start_pos - call_event.start_pos)
            end_dist = abs(truth_event.end_pos - call_event.end_pos)
            if start_dist > self.reference_distance or end_dist > self.reference_distance:
                return False

            # Check size similarity
            truth_size = abs(truth_event.end_pos - truth_event.start_pos)
            call_size = abs(call_event.end_pos - call_event.start_pos)
            if truth_size == 0 or call_size == 0:
                return False
            size_ratio = min(truth_size, call_size) / max(truth_size, call_size)
            if size_ratio < self.size_similarity:
                return False

            # Check sequence similarity if enabled
            if self.enable_sequence_comparison:
                similarity = self._calculate_sequence_similarity(truth_event, call_event)
                if similarity < self.sequence_similarity:
                    return False

            # Check reciprocal overlap
            overlap = self._calculate_overlap(truth_event, call_event)
            return not overlap < self.reciprocal_overlap

        except AttributeError as e:
            self.logger.warning(f"Error comparing events: {e!s}")
            return False

    def _compare_tra_events(self, truth_event, call_event) -> bool:
        """Special comparison logic for translocation events."""
        try:
            # Check chromosomes match
            if truth_event.start_chrom != call_event.start_chrom or truth_event.end_chrom != call_event.end_chrom:
                return False

            # Check positions within reference distance
            start_dist = abs(truth_event.start_pos - call_event.start_pos)
            end_dist = abs(truth_event.end_pos - call_event.end_pos)
            if start_dist > self.reference_distance or end_dist > self.reference_distance:
                return False

            # Check strand consistency if available
            if hasattr(truth_event, "strand") and hasattr(call_event, "strand"):
                if truth_event.strand != call_event.strand:
                    return False

            return True

        except AttributeError as e:
            self.logger.warning(f"Error comparing TRA events: {e!s}")
            return False

    def _calculate_sequence_similarity(self, truth_event, call_event) -> float:
        """Calculate sequence similarity between events."""
        truth_seq = self._get_sequence_from_event(truth_event)
        call_seq = self._get_sequence_from_event(call_event)

        # If no sequence information available, assume similarity
        if not truth_seq or not call_seq:
            return 1.0

        try:
            from Levenshtein import ratio

            return ratio(truth_seq, call_seq)
        except ImportError:
            self.logger.warning("Levenshtein package not available, using simple comparison")
            return float(truth_seq == call_seq)

    def _get_sequence_from_event(self, event) -> str | None:
        """Extract sequence information from an event."""
        try:
            if hasattr(event, "alt_seq") and event.alt_seq:
                return event.alt_seq

            if hasattr(event, "info"):
                seq = event.info.get("SVSEQ", "")
                if seq:
                    return seq
                seq = event.info.get("SEQ", "")
                if seq:
                    return seq

            if hasattr(event, "alt") and len(event.alt) > 1 and not event.alt.startswith("<"):
                return event.alt

            return None

        except AttributeError:
            return None

    def _calculate_overlap(self, event1, event2) -> float:
        """Calculate reciprocal overlap between two events."""
        try:
            if event1.sv_type == "TRA" or event2.sv_type == "TRA":
                return 1.0  # TRA events are compared by breakpoints only

            overlap_start = max(event1.start_pos, event2.start_pos)
            overlap_end = min(event1.end_pos, event2.end_pos)

            if overlap_start >= overlap_end:
                return 0.0

            overlap_length = overlap_end - overlap_start
            event1_length = event1.end_pos - event1.start_pos
            event2_length = event2.end_pos - event2.start_pos

            if event1_length == 0 or event2_length == 0:
                return 0.0

            overlap_ratio1 = overlap_length / event1_length
            overlap_ratio2 = overlap_length / event2_length

            return min(overlap_ratio1, overlap_ratio2)

        except AttributeError as e:
            self.logger.warning(f"Error calculating overlap: {e!s}")
            return 0.0

    def _compare_events(self):
        """Compare truth and call events to identify matches."""
        self.logger.info("Filtering events...")
        filtered_truth = self._filter_events(self.truth_events)
        filtered_calls = self._filter_events(self.call_events)

        self.logger.info("Comparing events...")
        tp_base, tp_call, fp = [], [], []
        matched_truth = set()
        matched_calls = set()

        # Compare each call against truth
        for call_event in filtered_calls:
            found_match = False
            for truth_event in filtered_truth:
                if truth_event in matched_truth:
                    continue

                if self._meets_matching_criteria(truth_event, call_event):
                    tp_call.append(call_event)
                    tp_base.append(truth_event)
                    matched_truth.add(truth_event)
                    matched_calls.add(call_event)
                    found_match = True
                    break

            if not found_match:
                fp.append(call_event)

        # Collect unmatched truth events as FN
        fn = [event for event in filtered_truth if event not in matched_truth]

        self.results = {"tp_base": tp_base, "tp_call": tp_call, "fp": fp, "fn": fn}

        self.logger.info(f"Found {len(tp_call)} true positives, {len(fp)} false positives, {len(fn)} false negatives")

    def _write_results(self):
        """Write benchmark results to output directory."""
        self.logger.info("Writing results...")
        self.output_dir.mkdir(parents=True, exist_ok=True)

        write_vcf(self.output_dir / "tp-base.vcf", self.results["tp_base"])
        write_vcf(self.output_dir / "tp-call.vcf", self.results["tp_call"])
        write_vcf(self.output_dir / "fp.vcf", self.results["fp"])
        write_vcf(self.output_dir / "fn.vcf", self.results["fn"])

        metrics = calculate_metrics(self.results)
        write_summary(self.output_dir / "summary.json", metrics)
