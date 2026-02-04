class EventTransformer:
    """Base class for transforming events."""

    def __init__(self, transform_strategies):
        self.transform_strategies = transform_strategies

    def apply_transforms(self, events):
        # The base implementation can be empty or provide a default behavior.
        pass

    def write_vcf(self, headers, transformed_events, output_file):
        # The base implementation can be empty or provide a default behavior.
        pass
