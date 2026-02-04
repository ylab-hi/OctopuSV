from .sv_selector import select_representative_sv
from .bnd_merge_logic import should_merge_bnd


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
        # Extract target chromosome from BND ALT field
        from .bnd_merge_logic import parse_bnd_alt
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

    def merge_events(self):
        """Merge nearly identical BND events for each chromosome pair.

        Returns:
            dict: A dictionary where:
                - Keys are chromosome pairs (tuple)
                - Values are lists of event groups (each group contains related events)

        The merging process:
        1. Iterates through each chromosome pair
        2. For each pair, groups nearly identical events based on:
           - Identical BND patterns (no reciprocal matching)
           - Breakpoint proximity (using delta)
           - Exact chromosome pair matching
        3. Maintains separate groups for non-matching events
        """
        all_chromosome_pair_events = {}
        for chromosome_pair, unmerged_events in self.bnd_events.items():
            merged_event_groups = []
            for current_event in unmerged_events:
                event_was_merged = False
                for _idx, event_group in enumerate(merged_event_groups):
                    existing_event = event_group[0]
                    if should_merge_bnd(existing_event, current_event, self.delta):
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
        all_chromosome_pair_events = self.merge_events()
        merged_events = []
        for _chromosome_pair, event_groups in all_chromosome_pair_events.items():
            for event_group in event_groups:
                representative_sv = select_representative_sv(event_group)
                # source_file merging is handled within select_representative_sv
                merged_events.append(representative_sv)
        return merged_events