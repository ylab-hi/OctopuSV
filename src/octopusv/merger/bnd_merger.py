from .sv_selector import select_representative_sv
from .bnd_merge_logic import parse_bnd_alt, should_merge_bnd


class BNDMerger:
    """A specialized merger class for handling breakend (BND) structural variants.

    Implements conservative clustering and merging of BND events between chromosome pairs.
    BND events require stricter matching criteria compared to TRA events, as they represent
    unclassified breakpoints that should be handled conservatively.

    The merger maintains a collection of BND events organized by chromosome pairs
    and provides methods to merge nearly identical breakends.
    """

    def __init__(self, delta=50):
        """Initialize the BND merger with configurable parameters.

        Args:
            delta (int): Distance threshold in base pairs for merging nearby breakpoints (default: 50)
                        This is more conservative than TRA merging to ensure precision

        The merger organizes BND events by chromosome pairs to efficiently handle
        breakpoints between specific chromosome combinations.
        """
        self.bnd_events = {}  # Dictionary to store BND events by chromosome pairs
        self.delta = delta
        self._merged_events_cache = None

    def add_event(self, event):
        """Add a BND event to the merger, organizing events by chromosome pairs.

        Args:
            event (SVCFEvent): A breakend event to be added to the merger

        The method:
        1. Extracts target chromosome from BND ALT field
        2. Creates a sorted tuple of chromosomes involved in the breakend
        3. Initializes a new list for previously unseen chromosome pairs
        4. Adds the event to the appropriate chromosome pair's event list
        """
        self._merged_events_cache = None

        # Extract target chromosome from BND ALT field
        _, target_chr, _ = parse_bnd_alt(event.alt)

        if target_chr is None:
            # If we can't parse the target chromosome, use the event's own chromosome
            # This handles cases where BND doesn't reference another chromosome
            key = (event.chrom, event.chrom)
        else:
            # Create sorted tuple for consistent key generation
            key = tuple(sorted([event.chrom, target_chr]))

        if key not in self.bnd_events:
            self.bnd_events[key] = []
        self.bnd_events[key].append(event)

    def _bnd_candidate_key(self, event):
        """Return the exact non-coordinate BND key required by should_merge_bnd().

        A valid key contains the BND pattern, local chromosome, and target
        chromosome. Events with incomplete ALT parsing return None and are kept
        on the original full-scan path.
        """
        pattern, target_chr, target_pos = parse_bnd_alt(event.alt)
        if target_chr is None or target_pos is None:
            return None, None
        return (pattern, event.chrom, target_chr), target_pos

    def merge_events(self):
        """Merge nearly identical BND events for each chromosome pair.

        Candidate pruning only removes comparisons that should_merge_bnd() is
        guaranteed to reject: mismatched BND pattern/chromosomes or breakpoint
        coordinates farther than delta. Final group assignment still uses the
        unchanged should_merge_bnd() function and original first-match order.

        Returns:
            dict: A dictionary where:
                - Keys are chromosome pairs (tuple)
                - Values are lists of event groups (each group contains related events)
        """
        all_chromosome_pair_events = {}

        for chromosome_pair, unmerged_events in self.bnd_events.items():
            merged_event_groups = []

            # Each group seed is indexed once. The index stores group indices,
            # never reordered groups, so candidate checks can be restored to the
            # exact original first-match order with sorted(indices).
            candidate_bins = {}

            # A non-negative delta is required for safe coordinate binning.
            # Negative values retain the original full-scan behavior.
            bin_width = self.delta + 1 if self.delta >= 0 else None

            for current_event in unmerged_events:
                event_was_merged = False
                candidate_key, target_pos = self._bnd_candidate_key(current_event)

                if candidate_key is None or bin_width is None:
                    candidate_indices = range(len(merged_event_groups))
                else:
                    local_bin = int(current_event.pos) // bin_width
                    target_bin = int(target_pos) // bin_width
                    candidate_indices = []

                    # If two coordinates differ by <= delta and bin width is
                    # delta + 1, their bins can differ by at most one.
                    for local_offset in (-1, 0, 1):
                        for target_offset in (-1, 0, 1):
                            bin_key = (
                                candidate_key,
                                local_bin + local_offset,
                                target_bin + target_offset,
                            )
                            candidate_indices.extend(candidate_bins.get(bin_key, ()))

                    candidate_indices.sort()

                for idx in candidate_indices:
                    existing_event = merged_event_groups[idx][0]
                    if should_merge_bnd(existing_event, current_event, self.delta):
                        merged_event_groups[idx].append(current_event)
                        event_was_merged = True
                        break

                if not event_was_merged:
                    group_index = len(merged_event_groups)
                    merged_event_groups.append([current_event])

                    # Only fully parsed seeds enter the pruning index. Incomplete
                    # records remain reachable through the full-scan fallback.
                    if candidate_key is not None and bin_width is not None:
                        local_bin = int(current_event.pos) // bin_width
                        target_bin = int(target_pos) // bin_width
                        bin_key = (candidate_key, local_bin, target_bin)
                        candidate_bins.setdefault(bin_key, []).append(group_index)

            all_chromosome_pair_events[chromosome_pair] = merged_event_groups

        return all_chromosome_pair_events

    def get_merged_events(self):
        """Get the final list of merged BND events with representative selection.

        Returns:
            list: A list of representative SVCFEvent objects, where each event:
                - Represents a group of nearly identical BND events
                - Is selected based on quality metrics (support, quality score, etc.)
                - Contains merged source file information

        Process:
        1. Obtains grouped events from merge_events()
        2. For each group, selects a representative event using select_representative_sv
        3. Ensures source file information is properly merged
        4. Returns the list of representative events
        """
        if self._merged_events_cache is not None:
            return self._merged_events_cache

        all_chromosome_pair_events = self.merge_events()
        merged_events = []
        for _chromosome_pair, event_groups in all_chromosome_pair_events.items():
            for event_group in event_groups:
                representative_sv = select_representative_sv(event_group)
                # source_file merging is handled within select_representative_sv
                merged_events.append(representative_sv)

        self._merged_events_cache = merged_events
        return self._merged_events_cache
