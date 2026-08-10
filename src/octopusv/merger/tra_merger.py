from .sv_selector import select_representative_sv
from .TRA_merge_logic import _breakend_orientation_signature, should_merge_tra


class TRAMerger:
    """A specialized merger class for handling translocation (TRA) structural variants.

    Implements clustering and merging of TRA events between chromosome pairs.

    The merger maintains a collection of TRA events organized by chromosome pairs
    and provides methods to merge overlapping translocations while considering
    strand consistency and overlap ratios.
    """

    def __init__(self, delta=50, min_overlap_ratio=0.5, strand_consistency=True):
        """Initialize the TRA merger with configurable parameters.

        Args:
            delta (int): Distance threshold in base pairs for merging nearby breakpoints (default: 50)
            min_overlap_ratio (float): Minimum required overlap ratio between events to be merged (default: 0.5)
            strand_consistency (bool): Whether to enforce strand consistency during merging (default: True)

        The merger organizes TRA events by chromosome pairs to efficiently handle
        translocations between specific chromosome combinations.
        """
        self.tra_events = {}  # Dictionary to store TRA events by chromosome pairs
        self.delta = delta
        self.min_overlap_ratio = min_overlap_ratio
        self.strand_consistency = strand_consistency
        self._merged_events_cache = None

    def add_event(self, event):
        """Add a TRA event to the merger, organizing events by chromosome pairs.

        Args:
            event (SVCFEvent): A translocation event to be added to the merger

        The method:
        1. Creates a sorted tuple of chromosomes involved in the translocation
        2. Initializes a new list for previously unseen chromosome pairs
        3. Adds the event to the appropriate chromosome pair's event list
        """
        self._merged_events_cache = None

        key = tuple(sorted([event.start_chrom, event.end_chrom]))
        if key not in self.tra_events:
            self.tra_events[key] = []
        self.tra_events[key].append(event)

    @staticmethod
    def _adjacency_candidate_info(event):
        """Return a safe physical-adjacency candidate key and aligned positions.

        The key contains only requirements already enforced by
        should_merge_tra(): reverse-complement state plus the unordered pair of
        chromosome-side endpoints. When the two endpoint identities are
        distinct, their positions can also be aligned canonically for safe
        coordinate binning. Duplicate endpoint identities are left unbinned so
        reciprocal pairing remains fully delegated to should_merge_tra().
        """
        signature = _breakend_orientation_signature(event)
        if signature is None:
            return None, None

        local_endpoint, remote_endpoint, reverse_complement = signature
        local_identity = (local_endpoint[0], local_endpoint[2])
        remote_identity = (remote_endpoint[0], remote_endpoint[2])
        endpoint_sides = tuple(sorted((local_identity, remote_identity)))
        key = (endpoint_sides, reverse_complement)

        if local_identity == remote_identity:
            return key, None

        if local_identity < remote_identity:
            positions = (int(local_endpoint[1]), int(remote_endpoint[1]))
        else:
            positions = (int(remote_endpoint[1]), int(local_endpoint[1]))

        return key, positions

    def merge_events(self):
        """Merge overlapping TRA events for each chromosome pair.

        Candidate pruning only removes comparisons that the unchanged
        should_merge_tra() function is guaranteed to reject.

        With strand consistency enabled, parseable physical-adjacency signatures
        are first separated by the same chromosome-side and reverse-complement
        requirements used by should_merge_tra(). For distinct endpoint
        identities, candidates are additionally limited to the same two-breakpoint
        tolerance window. Signature-less groups remain candidates because the
        original logic falls back to legacy coordinate matching whenever either
        signature is unavailable.

        Candidate group indices are checked in original creation order, so the
        existing first-match behavior is preserved exactly.
        """
        all_chromosome_pair_events = {}

        for chromosome_pair, unmerged_events in self.tra_events.items():
            merged_event_groups = []

            adjacency_group_indices = {}
            adjacency_bins = {}
            fallback_group_indices = []

            tra_delta = self.delta * 2
            bin_width = tra_delta + 1 if tra_delta >= 0 else None

            for current_event in unmerged_events:
                event_was_merged = False

                if not self.strand_consistency or bin_width is None:
                    current_key = None
                    current_positions = None
                    candidate_indices = range(len(merged_event_groups))
                else:
                    current_key, current_positions = self._adjacency_candidate_info(current_event)

                    if current_key is None:
                        # should_merge_tra() uses legacy coordinate matching when
                        # either event lacks a valid signature, so no signature-
                        # based pruning is safe for this current event.
                        candidate_indices = range(len(merged_event_groups))
                    elif current_positions is None:
                        # Duplicate chromosome-side endpoint identities make the
                        # same/reciprocal positional pairing ambiguous. Preserve
                        # all groups with the same physical-adjacency key.
                        candidate_indices = sorted(
                            adjacency_group_indices.get(current_key, [])
                            + fallback_group_indices
                        )
                    else:
                        pos1_bin = current_positions[0] // bin_width
                        pos2_bin = current_positions[1] // bin_width
                        candidate_indices = list(fallback_group_indices)

                        # If both aligned breakpoint differences are <= tra_delta
                        # and bin width is tra_delta + 1, each bin can differ by
                        # at most one from the matching group seed.
                        for pos1_offset in (-1, 0, 1):
                            for pos2_offset in (-1, 0, 1):
                                bin_key = (
                                    current_key,
                                    pos1_bin + pos1_offset,
                                    pos2_bin + pos2_offset,
                                )
                                candidate_indices.extend(adjacency_bins.get(bin_key, ()))

                        candidate_indices.sort()

                for idx in candidate_indices:
                    existing_event = merged_event_groups[idx][0]
                    if should_merge_tra(
                        existing_event,
                        current_event,
                        self.delta,
                        self.min_overlap_ratio,
                        self.strand_consistency,
                    ):
                        merged_event_groups[idx].append(current_event)
                        event_was_merged = True
                        break

                if not event_was_merged:
                    group_index = len(merged_event_groups)
                    merged_event_groups.append([current_event])

                    if self.strand_consistency and bin_width is not None:
                        if current_key is None:
                            fallback_group_indices.append(group_index)
                        else:
                            adjacency_group_indices.setdefault(current_key, []).append(group_index)
                            if current_positions is not None:
                                pos1_bin = current_positions[0] // bin_width
                                pos2_bin = current_positions[1] // bin_width
                                bin_key = (current_key, pos1_bin, pos2_bin)
                                adjacency_bins.setdefault(bin_key, []).append(group_index)

            all_chromosome_pair_events[chromosome_pair] = merged_event_groups

        return all_chromosome_pair_events

    def get_merged_events(self):
        """Get the final list of merged TRA events with representative selection.

        Returns:
            list: A list of representative SVCFEvent objects, where each event:
                - Represents a group of overlapping TRA events
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


"""
self.tra_events = {
    ('chr1', 'chr2'): [
        (1000, 2000, 'FileA'),
        (1050, 2030, 'FileB'),
        (1500, 2500, 'FileC'),
        (1520, 2510, 'FileD')
    ],
    ('chr3', 'chr4'): [
        (3000, 4000, 'FileA'),
        (3100, 4100, 'FileB')
    ]
}

self.distance_threshold = 100

all_chromosome_pair_events = {
    ('chr1', 'chr2'): [
        [1000, 2000, {'FileA', 'FileB'}],
        [1500, 2500, {'FileC', 'FileD'}]
    ],
    ('chr3', 'chr4'): [
        [3000, 4000, {'FileA', 'FileB'}]
    ]
}

merged_events_for_certain_chrom_pair = [
    [1000, 2000, {'FileA', 'FileB'}],
    [1500, 2500, {'FileC', 'FileD'}]
]

"""
