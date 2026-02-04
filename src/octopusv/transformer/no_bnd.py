from .base import EventTransformer


class NonBNDTransformer(EventTransformer):
    """Transformer for non-BND events. Currently, the transformation is a no-op.

    but the structure is in place for future transformation strategies.
    """

    def apply_transforms(self, events):
        """Apply all transformation strategies to a list of events, even if they do nothing."""
        for event in events:
            for strategy in self.transform_strategies:
                # Each strategy is applied to each event, regardless of whether it transforms the event.
                strategy.convert(event)
        return events
