from .bnd_merger import BNDMerger
from .sv_merge_logic import should_merge
from .sv_merge_selection import MergeSelectionMixin
from .sv_merge_writer import MergeWriterMixin
from .sv_selector import select_representative_sv
from .tra_merger import TRAMerger


class SVMerger(MergeSelectionMixin, MergeWriterMixin):
    def __init__(
            self,
            classified_events,
            all_input_files,
            tra_delta=50,
            tra_min_overlap_ratio=0.5,
            tra_strand_consistency=True,
            max_distance=50,
            max_length_ratio=1.3,
            min_jaccard=0.7,
            bnd_delta=50,
    ):
        """Initialize SVMerger with the given parameters and events.

        Args:
            classified_events: Dictionary of classified SV events
            all_input_files: List of all input file paths
            tra_delta: Position uncertainty threshold for TRA events
            tra_min_overlap_ratio: Minimum overlap ratio for TRA events
            tra_strand_consistency: Whether to require strand consistency for TRA events
            max_distance: Maximum allowed distance between start or end positions
            max_length_ratio: Maximum allowed ratio between event lengths
            min_jaccard: Minimum required Jaccard index for overlap
            bnd_delta: Position uncertainty threshold for BND events (default: 50)
        """
        self.classified_events = classified_events
        self.all_input_files = [str(file) for file in all_input_files]
        self.merged_events: dict[str, dict[str, list]] = {}
        self.event_groups: dict[str, dict[str, list[list]]] = {}
        self.tra_merger = TRAMerger(tra_delta, tra_min_overlap_ratio, tra_strand_consistency)
        self.bnd_merger = BNDMerger(bnd_delta)
        self.max_distance = max_distance
        self.max_length_ratio = max_length_ratio
        self.min_jaccard = min_jaccard

    def merge(self):
        """Merge all SV events based on their types and chromosomes."""
        for sv_type, chromosomes in self.classified_events.items():
            if sv_type == "TRA":
                for (_chr1, _chr2), events in chromosomes.items():
                    for event in events:
                        self.tra_merger.add_event(event)
            elif sv_type == "BND":
                for chromosome, events in chromosomes.items():
                    for event in events:
                        self.bnd_merger.add_event(event)
            else:
                if sv_type not in self.merged_events:
                    self.merged_events[sv_type] = {}
                    self.event_groups[sv_type] = {}
                for chromosome, events in chromosomes.items():
                    if chromosome not in self.merged_events[sv_type]:
                        self.merged_events[sv_type][chromosome] = []
                        self.event_groups[sv_type][chromosome] = []
                    for event in sorted(events, key=lambda e: (e.start_pos, e.end_pos, e.sv_id)):
                        self.add_and_merge_event(sv_type, chromosome, event)

    def add_and_merge_event(self, sv_type, chromosome, new_event):
        """Add a new event and merge it with existing events if possible."""
        events = self.merged_events[sv_type][chromosome]
        event_groups = self.event_groups[sv_type][chromosome]
        for idx, existing_event in enumerate(events):
            if should_merge(existing_event, new_event, self.max_distance, self.max_length_ratio, self.min_jaccard):
                event_groups[idx].append(new_event)
                return
        events.append(new_event)
        event_groups.append([new_event])

    def get_events(self, sv_type, chromosome, start, end):
        """Get events of given type within the specified region."""
        if sv_type == "TRA":
            return self.tra_merger.get_merged_events()
        elif sv_type == "BND":
            return self.bnd_merger.get_merged_events()
        if sv_type in self.event_groups and chromosome in self.event_groups[sv_type]:
            events = []
            for sv_group in self.event_groups[sv_type][chromosome]:
                representative_sv = select_representative_sv(sv_group)
                if representative_sv.start_pos <= end and representative_sv.end_pos >= start:
                    events.append(representative_sv)
            return events
        return []

    def get_all_merged_events(self):
        """Get all merged events across all types and chromosomes."""
        merged_events = []
        merged_events.extend(self.tra_merger.get_merged_events())
        merged_events.extend(self.bnd_merger.get_merged_events())
        for sv_type in self.event_groups:
            for chromosome in self.event_groups[sv_type]:
                for sv_group in self.event_groups[sv_type][chromosome]:
                    representative_sv = select_representative_sv(sv_group)
                    merged_events.append(representative_sv)
        return merged_events
