class SVClassifiedByChromosome:
    """Further classifies SV events by chromosome within each SV type category.

    This helps reduce the search and processing space, especially when performing interval-based operations.

    Attributes:
        events_by_type (dict): A dictionary with SV types as keys and lists of events as values.

    Output will be something like:
    {
    "INS": {
        "chr1": [event1, event2],
        "chr5": [event3]
    },
    "DEL": {
        "chr5": [event4, event5],
        "chr8": [event6]
    },
    "TRA": {
        ("chr10", "chr12"): [event7],
        ("chr1", "chr10"): [event8, event9]
    },
    "DUP": {
        "chr1": [event10, event11, event12],
        "chr7": [event13]
    },
    "INV": {
        "chr10": [event14],
        "chr12": [event15, event16],
        "chr5": [event17, event18]
    }
    }

    """

    def __init__(self, events_by_type):
        """Constructor to initialize the classifier.

        Parameters:
            events_by_type (dict): A dictionary containing events classified by SV type.
        """
        self.events_by_type = events_by_type
        # Initialize an empty dictionary for each SV type to store events classified by chromosome.
        self.classified_by_chromosome = {
            sv_type: {} for sv_type in events_by_type
        }  # This output is a double layer dictionary.

    def classify(self):
        """Classifies events within each SV type by chromosome, handling TRA differently by considering both start and end chromosomes."""
        for sv_type, events in self.events_by_type.items():
            if sv_type == "TRA":
                # Special handling for TRA events, classified by both start and end chromosomes.
                self._classify_tra(events)
            else:
                # Classify regular SV types by the start chromosome.
                self._classify_regular(events, sv_type)

    def _classify_regular(self, events, sv_type):
        """Classifies regular SV types by start chromosome.

        Parameters:
            events (list): List of events.
            sv_type (str): Current SV type being processed.
        """
        for event in events:
            chrom_key = event.start_chrom
            # Ensure each chromosome key is in the dictionary, initialize an empty list if not.
            if chrom_key not in self.classified_by_chromosome[sv_type]:
                self.classified_by_chromosome[sv_type][chrom_key] = []
            # Append the event under its corresponding chromosome key.
            self.classified_by_chromosome[sv_type][chrom_key].append(event)

    def _classify_tra(self, events):
        """Specially handles TRA events by classifying them using both start and end chromosomes.

        Parameters:
            events (list): List of TRA events.
        """
        for event in events:
            # Use a tuple of start and end chromosomes as a key.
            chrom_key = (event.start_chrom, event.end_chrom)
            # Ensure the chromosome key exists in the dictionary, initialize an empty list if not.
            if chrom_key not in self.classified_by_chromosome["TRA"]:
                self.classified_by_chromosome["TRA"][chrom_key] = []
            # Append the event under its corresponding chromosome pair key.
            self.classified_by_chromosome["TRA"][chrom_key].append(event)

    def get_classified_events(self):
        """Returns the dictionary containing SV events classified by chromosome.

        Returns:
            dict: A nested dictionary with SV types as top-level keys, chromosome pairs or single chromosomes as second-level keys,
            and lists of events as values.
        """
        return self.classified_by_chromosome
