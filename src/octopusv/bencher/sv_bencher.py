import logging
from pathlib import Path

from octopusv.utils.svcf_parser import SVCFFileEventCreator
from octopusv.utils.svcf_validator import validate_versioned_svcf_for_downstream

from .bench_utils import calculate_metrics, read_safe_vcf_meta, write_summary, write_vcf


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
        """Validate versioned inputs, then parse truth and call events."""
        validate_versioned_svcf_for_downstream(
            self.truth_file, consumer="octopusv benchmark truth input"
        )
        validate_versioned_svcf_for_downstream(
            self.call_file, consumer="octopusv benchmark call input"
        )

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

            # Breakpoint events (TRA/BND) are defined by two genomic
            # breakpoints, potentially on different chromosomes. Numeric
            # subtraction between their coordinates is not an SV length and
            # must never drive --size-min/--size-max filtering.
            try:
                if not self._is_breakpoint_event(event):
                    size = abs(event.end_pos - event.start_pos)
                    if size < self.size_min or size > self.size_max:
                        continue

                filtered.append(event)
            except AttributeError as e:
                self.logger.warning(f"Skipping malformed event: {e!s}")
                continue

        return filtered

    def _meets_matching_criteria(self, truth_event, call_event) -> bool:
        """Check if two events meet the matching criteria based on GIAB standards."""
        try:
            # Check SV type unless ignored. --type-ignore relaxes the exact
            # SVTYPE label, but it must not mix incompatible geometries: a
            # two-breakpoint event (TRA/BND) is never an interval event.
            if not self.type_ignore and truth_event.sv_type != call_event.sv_type:
                return False

            truth_is_breakpoint = self._is_breakpoint_event(truth_event)
            call_is_breakpoint = self._is_breakpoint_event(call_event)
            if truth_is_breakpoint != call_is_breakpoint:
                return False

            if truth_is_breakpoint:
                return self._compare_breakpoint_events(truth_event, call_event)

            # Ordinary interval SVs must occur on the same chromosome.
            # end_chrom is intentionally not required here because legacy
            # SVCF may carry CHR2=. for otherwise valid intra-chromosomal
            # records; CHROM/POS and END define the interval geometry.
            if truth_event.start_chrom != call_event.start_chrom:
                return False

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

    @staticmethod
    def _is_breakpoint_event(event) -> bool:
        """Return whether an event uses two-breakpoint TRA/BND geometry."""
        return event.sv_type in {"TRA", "BND"}

    def _compare_breakpoint_events(self, truth_event, call_event) -> bool:
        """Compare TRA/BND events by their ordered breakpoint geometry.

        Endpoint swapping is intentionally not attempted here; this preserves
        the historical benchmark contract while fixing chromosome-aware BND
        handling.
        """
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
            self.logger.warning(f"Error comparing breakpoint events: {e!s}")
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
            if self._is_breakpoint_event(event1) or self._is_breakpoint_event(event2):
                return 1.0  # TRA/BND events are compared by breakpoints only

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

        truth_meta = read_safe_vcf_meta(self.truth_file)
        call_meta = read_safe_vcf_meta(self.call_file)

        write_vcf(
            self.output_dir / "tp-base.vcf",
            self.results["tp_base"],
            source_meta_lines=truth_meta,
        )
        write_vcf(
            self.output_dir / "tp-call.vcf",
            self.results["tp_call"],
            source_meta_lines=call_meta,
        )
        write_vcf(
            self.output_dir / "fp.vcf",
            self.results["fp"],
            source_meta_lines=call_meta,
        )
        write_vcf(
            self.output_dir / "fn.vcf",
            self.results["fn"],
            source_meta_lines=truth_meta,
        )

        metrics = calculate_metrics(self.results)
        write_summary(self.output_dir / "summary.json", metrics)
