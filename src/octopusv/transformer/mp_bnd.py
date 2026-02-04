from .base import EventTransformer


class MatePairBNDTransformer(EventTransformer):
    """Class for transforming mate pair BND events."""

    def apply_transforms(self, mate_pairs):
        transformed_events = []
        for pair in mate_pairs:  # Pair is a tuple of mate events.
            for strategy in self.transform_strategies:
                transformed_events.extend(
                    strategy.convert(pair),
                )  # Strategy must be able to handle a tuple of events.
        return transformed_events
