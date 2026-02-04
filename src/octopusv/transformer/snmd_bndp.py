from .base import EventTransformer


class SpecialNoMateDiffBNDPairTransformer(EventTransformer):
    """Class for transforming special_no_mate_diff_bnd_pair events."""

    def apply_transforms(self, event_pairs):
        transformed_events = []
        for pair in event_pairs:  # Pair is a tuple of two events.
            for strategy in self.transform_strategies:
                transformed_events.extend(
                    strategy.convert(pair),
                )  # Strategy must be able to handle a tuple of events.
        return transformed_events
