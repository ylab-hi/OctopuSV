from .sv_selector import select_representative_sv
from .TRA_merge_logic import should_merge_tra


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

    def add_event(self, event):
        """Add a TRA event to the merger, organizing events by chromosome pairs.

        Args:
            event (SVCFEvent): A translocation event to be added to the merger

        The method:
        1. Creates a sorted tuple of chromosomes involved in the translocation
        2. Initializes a new list for previously unseen chromosome pairs
        3. Adds the event to the appropriate chromosome pair's event list
        """
        key = tuple(sorted([event.start_chrom, event.end_chrom]))
        if key not in self.tra_events:
            self.tra_events[key] = []
        self.tra_events[key].append(event)

    def merge_events(self):
        """Merge overlapping TRA events for each chromosome pair.

        Returns:
            dict: A dictionary where:
                - Keys are chromosome pairs (tuple)
                - Values are lists of event groups (each group contains related events)

        The merging process:
        1. Iterates through each chromosome pair
        2. For each pair, groups overlapping events based on:
           - Breakpoint proximity (using delta)
           - Overlap ratio threshold
           - Strand consistency (if enabled)
        3. Maintains separate groups for non-overlapping events
        """
        all_chromosome_pair_events = {}
        for chromosome_pair, unmerged_events in self.tra_events.items():
            merged_event_groups = []
            for current_event in unmerged_events:
                event_was_merged = False
                for _idx, event_group in enumerate(merged_event_groups):
                    existing_event = event_group[0]
                    if should_merge_tra(
                        existing_event, current_event, self.delta, self.min_overlap_ratio, self.strand_consistency
                    ):
                        # Add current event to existing group if it meets merge criteria
                        event_group.append(current_event)
                        event_was_merged = True
                        break
                if not event_was_merged:
                    # Create new group if event doesn't match any existing groups
                    merged_event_groups.append([current_event])
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
        all_chromosome_pair_events = self.merge_events()
        merged_events = []
        for _chromosome_pair, event_groups in all_chromosome_pair_events.items():
            for event_group in event_groups:
                representative_sv = select_representative_sv(event_group)
                # source_file merging is handled within select_representative_sv
                merged_events.append(representative_sv)
        return merged_events


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
