from .base import EventTransformer


class MatePairBNDTransformer(EventTransformer):
    """Class for transforming mate pair BND events.

    The strategy tables are mutually exclusive: every mate pair is converted by
    exactly one strategy.  A pair matched by more than one strategy is an
    internal error (it used to be emitted twice).  A pair matched by none is
    kept in ``self.unmatched`` so the caller can route both records through
    the single-TRA path instead of dropping them.
    """

    def __init__(self, transform_strategies):
        super().__init__(transform_strategies)
        self.unmatched = []
        self.same_adjacency_pairs = 0
        self.distinct_adjacency_pairs = 0

    def apply_transforms(self, mate_pairs):
        transformed_events = []
        self.unmatched = []
        for pair in mate_pairs:  # Pair is a tuple of mate events.
            results = [(strategy, strategy.convert(pair)) for strategy in self.transform_strategies]
            matched = [(strategy, events) for strategy, events in results if events]
            if len(matched) > 1:
                names = ", ".join(type(strategy).__name__ for strategy, _ in matched)
                raise RuntimeError(
                    f"BND mate pair {pair[0].id} / {pair[1].id} matched several conversion rules ({names})"
                )
            if not matched:
                self.unmatched.append(pair)
                continue
            events = matched[0][1]
            if len(events) == 1:
                self.same_adjacency_pairs += 1
            else:
                self.distinct_adjacency_pairs += 1
            transformed_events.extend(events)
        return transformed_events
