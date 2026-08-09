import os
import re


class MergeSelectionMixin:
    def get_events_by_source(self, sources, operation="union"):
        """Get events based on their source files and specified operation."""
        tra_events = self.tra_merger.get_merged_events()
        bnd_events = self.bnd_merger.get_merged_events()
        other_events = self.get_all_merged_events()
        other_events = [e for e in other_events if e.sv_type not in ["TRA", "BND"]]
        sources_set = {os.path.basename(source) for source in sources}

        if operation == "union":
            tra_filtered = [
                event for event in tra_events
                if sources_set.intersection({os.path.basename(s) for s in event.source_file.split(",")})
            ]
            bnd_filtered = [
                event for event in bnd_events
                if sources_set.intersection({os.path.basename(s) for s in event.source_file.split(",")})
            ]
            other_filtered = [
                event for event in other_events
                if sources_set.intersection({os.path.basename(s) for s in event.source_file.split(",")})
            ]
        elif operation == "intersection":
            tra_filtered = [
                event for event in tra_events
                if sources_set.issubset({os.path.basename(s) for s in event.source_file.split(",")})
            ]
            bnd_filtered = [
                event for event in bnd_events
                if sources_set.issubset({os.path.basename(s) for s in event.source_file.split(",")})
            ]
            other_filtered = [
                event for event in other_events
                if sources_set.issubset({os.path.basename(s) for s in event.source_file.split(",")})
            ]
        elif operation == "specific":
            source_file = next(iter(sources_set))
            other_files = {os.path.basename(f) for f in self.all_input_files} - {source_file}
            tra_filtered = [
                event for event in tra_events
                if source_file in [os.path.basename(s) for s in event.source_file.split(",")] and not any(
                    other in [os.path.basename(s) for s in event.source_file.split(",")] for other in other_files
                )
            ]
            bnd_filtered = [
                event for event in bnd_events
                if source_file in [os.path.basename(s) for s in event.source_file.split(",")] and not any(
                    other in [os.path.basename(s) for s in event.source_file.split(",")] for other in other_files
                )
            ]
            other_filtered = [
                event for event in other_events
                if source_file in [os.path.basename(s) for s in event.source_file.split(",")] and not any(
                    other in [os.path.basename(s) for s in event.source_file.split(",")] for other in other_files
                )
            ]
        else:
            msg = f"Unsupported operation: {operation}"
            raise ValueError(msg)

        return other_filtered + tra_filtered + bnd_filtered

    def get_events_by_exact_support(self, exact_support):
        """Get events supported by exactly N files."""
        tra_events = self.tra_merger.get_merged_events()
        bnd_events = self.bnd_merger.get_merged_events()
        other_events = self.get_all_merged_events()
        other_events = [e for e in other_events if e.sv_type not in ["TRA", "BND"]]

        tra_filtered = [event for event in tra_events if len(set(event.source_file.split(","))) == exact_support]
        bnd_filtered = [event for event in bnd_events if len(set(event.source_file.split(","))) == exact_support]
        other_filtered = [event for event in other_events if len(set(event.source_file.split(","))) == exact_support]
        return other_filtered + tra_filtered + bnd_filtered

    def get_events_by_support_range(self, min_support=None, max_support=None):
        """Get events supported by a range of files."""
        tra_events = self.tra_merger.get_merged_events()
        bnd_events = self.bnd_merger.get_merged_events()
        other_events = self.get_all_merged_events()
        other_events = [e for e in other_events if e.sv_type not in ["TRA", "BND"]]

        def within_range(event):
            support_count = len(set(event.source_file.split(",")))
            if min_support is not None and support_count < min_support:
                return False
            return not (max_support is not None and support_count > max_support)

        tra_filtered = [event for event in tra_events if within_range(event)]
        bnd_filtered = [event for event in bnd_events if within_range(event)]
        other_filtered = [event for event in other_events if within_range(event)]
        return other_filtered + tra_filtered + bnd_filtered

    def get_events_by_expression(self, expression):
        """Get events that satisfy a logical expression."""
        tra_events = self.tra_merger.get_merged_events()
        bnd_events = self.bnd_merger.get_merged_events()
        other_events = self.get_all_merged_events()
        other_events = [e for e in other_events if e.sv_type not in ["TRA", "BND"]]

        tra_filtered = [
            event for event in tra_events
            if self.evaluate_expression(expression, [os.path.basename(s) for s in event.source_file.split(",")])
        ]
        bnd_filtered = [
            event for event in bnd_events
            if self.evaluate_expression(expression, [os.path.basename(s) for s in event.source_file.split(",")])
        ]
        other_filtered = [
            event for event in other_events
            if self.evaluate_expression(expression, [os.path.basename(s) for s in event.source_file.split(",")])
        ]
        return other_filtered + tra_filtered + bnd_filtered

    def evaluate_expression(self, expression, event_sources):
        """Evaluate a logical expression against event sources."""

        def make_identifier(file_path):
            file_name = os.path.basename(file_path)
            return re.sub(r"\W|^(?=\d)", "_", file_name)

        context = {}
        all_sources = set(self.all_input_files)
        for source in all_sources:
            identifier = make_identifier(source)
            context[identifier] = False

        for source in event_sources:
            identifier = make_identifier(source)
            context[identifier] = True

        expr = expression
        for source in all_sources:
            identifier = make_identifier(source)
            expr = re.sub(r"\b" + re.escape(os.path.basename(source)) + r"\b", identifier, expr)

        expr = expr.replace("AND", "and").replace("OR", "or").replace("NOT", "not")

        try:
            return eval(expr, {"__builtins__": {}}, context)
        except Exception as e:
            msg = f"Invalid expression: {e}"
            raise ValueError(msg)
