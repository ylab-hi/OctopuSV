class SVClassifierByType:
    """Classifies SVCF events by their SV type into separate categories.

    Attributes:
        events (list): List of SVCFEvent objects to be classified.
        classified_events (dict): Dictionary where keys are SV types and values are lists of events.
    """

    def __init__(self, events):
        self.events = events
        # Predefine categories with empty lists - now includes BND
        self.classified_events = {"INS": [], "DEL": [], "DUP": [], "INV": [], "TRA": [], "BND": []}
        self.valid_types = {"INS", "DEL", "DUP", "INV", "TRA", "BND"}  # Added BND to valid types

    def classify(self):
        """Classifies events based on their SV type. Only includes events whose types are in the valid_types set.

        Raises an error if an event has an undefined or invalid SV type.
        """
        for event in self.events:
            if event.sv_type not in self.valid_types:
                msg = f"Invalid SV type encountered: {event.sv_type}"
                raise ValueError(msg)

            # Directly append to the predefined list for each SV type
            self.classified_events[event.sv_type].append(event)

    def get_classified_events(self):
        """Returns the dictionary containing classified events.

        Returns:
            dict: A dictionary with SV types as keys and lists of corresponding events as values.
        """
        return self.classified_events