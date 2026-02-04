from .base import EventTransformer


class SingleTRATransformer(EventTransformer):
    """Class for transforming other_single_TRA."""

    def apply_transforms(self, events):
        for event in events:  # Apply the conversion for each event.
            for strategy in self.transform_strategies:
                strategy.convert(event)
        return events  # Return the original list of events, which have been modified in-place.
