import ast
import os
import re

from octopusv.utils.source_path import normalize_source_path


def _expression_identifier(file_path):
    file_name = os.path.basename(str(file_path))
    return re.sub(r"\W|^(?=\d)", "_", file_name)


def _build_expression_source_entries(input_files, normalize_func=normalize_source_path):
    """Build unambiguous basename/identifier mappings for --expression."""
    by_basename = {}
    by_identifier = {}
    entries = []

    for source in input_files:
        source_text = str(source)
        normalized = normalize_func(source_text)
        basename = os.path.basename(source_text)
        identifier = _expression_identifier(basename)

        previous = by_basename.get(basename)
        if previous is not None and previous != normalized:
            raise ValueError(
                "Cannot use --expression because multiple merge inputs "
                f"share the basename {basename!r}. Rename the files so "
                "each expression source name is unique."
            )
        by_basename[basename] = normalized

        previous_basename = by_identifier.get(identifier)
        if previous_basename is not None and previous_basename != basename:
            raise ValueError(
                "Cannot use --expression because input basenames "
                f"{previous_basename!r} and {basename!r} map to the same "
                f"expression identifier {identifier!r}. Rename one file."
            )
        by_identifier[identifier] = basename

        entries.append((basename, identifier, normalized))

    return entries


def _compile_expression(expression, entries):
    """Replace source basenames once and compile one logical expression."""
    expr = str(expression)

    # Replace longer basenames first in case one filename is a prefix of
    # another. The public expression syntax remains basename-based.
    for basename, identifier, _normalized in sorted(
        entries,
        key=lambda entry: len(entry[0]),
        reverse=True,
    ):
        expr = re.sub(
            r"(?<!\w)" + re.escape(basename) + r"(?!\w)",
            identifier,
            expr,
        )

    expr = re.sub(r"\bAND\b", "and", expr)
    expr = re.sub(r"\bOR\b", "or", expr)
    expr = re.sub(r"\bNOT\b", "not", expr)

    try:
        tree = ast.parse(expr, mode="eval")
    except SyntaxError as exc:
        detail = exc.msg or str(exc)
        raise ValueError(f"Invalid expression: {detail}") from exc

    allowed_identifiers = {identifier for _basename, identifier, _normalized in entries}
    referenced_identifiers = {
        node.id
        for node in ast.walk(tree)
        if isinstance(node, ast.Name)
    }
    unknown_identifiers = referenced_identifiers - allowed_identifiers
    if unknown_identifiers:
        unknown = ", ".join(sorted(unknown_identifiers))
        raise ValueError(
            "Invalid expression: unknown source identifier(s): "
            f"{unknown}. Use input file basenames in --expression."
        )

    return compile(tree, "<octopusv-expression>", "eval"), referenced_identifiers


def validate_selection_inputs(*, input_files, specific=None, expression=None):
    """Validate selection arguments before expensive parsing/merging starts."""
    if specific:
        normalized_inputs = {
            normalize_source_path(str(source))
            for source in input_files
        }
        for source in specific:
            source_text = str(source)
            normalized = normalize_source_path(source_text)
            if normalized not in normalized_inputs:
                raise ValueError(
                    "Requested source is not one of the merge inputs: "
                    f"{source_text!r}."
                )

    if expression:
        entries = _build_expression_source_entries(input_files)
        _compile_expression(expression, entries)


class MergeSelectionMixin:
    def _normalize_source(self, source):
        """Normalize one source identity once per distinct raw source string."""
        source_text = str(source)
        cache = getattr(self, "_source_identity_cache", None)
        if cache is None:
            cache = {}
            self._source_identity_cache = cache

        if source_text not in cache:
            cache[source_text] = normalize_source_path(source_text)
        return cache[source_text]

    def _normalized_input_sources(self):
        """Return exact normalized identities for all merge inputs."""
        normalized = {}
        for source in self.all_input_files:
            source_text = str(source)
            key = self._normalize_source(source_text)
            previous = normalized.get(key)
            if previous is not None and previous != source_text:
                raise ValueError(
                    "Merge inputs resolve to the same physical file: "
                    f"{previous!r} and {source_text!r}."
                )
            normalized[key] = source_text
        return normalized

    def _normalize_requested_sources(self, sources):
        """Normalize requested sources and ensure they belong to this merge."""
        input_sources = self._normalized_input_sources()
        requested = set()

        for source in sources:
            source_text = str(source)
            normalized = self._normalize_source(source_text)
            if normalized not in input_sources:
                raise ValueError(
                    "Requested source is not one of the merge inputs: "
                    f"{source_text!r}."
                )
            requested.add(normalized)

        return requested

    def _event_source_paths(self, event):
        """Return exact normalized source identities carried by one event."""
        return {
            self._normalize_source(source.strip())
            for source in str(getattr(event, "source_file", "")).split(",")
            if source.strip()
        }

    def _all_selectable_events(self):
        """Return ordinary, TRA, and BND events in the historical output order."""
        tra_events = self.tra_merger.get_merged_events()
        bnd_events = self.bnd_merger.get_merged_events()
        other_events = self.get_all_merged_events()
        other_events = [e for e in other_events if e.sv_type not in ["TRA", "BND"]]
        return other_events, tra_events, bnd_events

    def get_events_by_source(self, sources, operation="union"):
        """Get events based on exact input-source identity.

        ``sources`` must resolve to merge inputs. Source identity is based on
        normalized real paths, never basenames, so same-named files from
        different directories remain distinct.
        """
        requested_sources = self._normalize_requested_sources(sources)
        other_events, tra_events, bnd_events = self._all_selectable_events()

        def selected(event):
            event_sources = self._event_source_paths(event)

            if operation == "union":
                return bool(requested_sources.intersection(event_sources))

            if operation == "intersection":
                return requested_sources.issubset(event_sources)

            if operation == "specific":
                # "Specifically supported by provided files" means every
                # supporting source belongs to the requested set, while at
                # least one requested source actually supports the event.
                return bool(event_sources) and event_sources.issubset(requested_sources)

            msg = f"Unsupported operation: {operation}"
            raise ValueError(msg)

        other_filtered = [event for event in other_events if selected(event)]
        tra_filtered = [event for event in tra_events if selected(event)]
        bnd_filtered = [event for event in bnd_events if selected(event)]
        return other_filtered + tra_filtered + bnd_filtered

    def get_events_by_exact_support(self, exact_support):
        """Get events supported by exactly N distinct physical input sources."""
        other_events, tra_events, bnd_events = self._all_selectable_events()

        def has_exact_support(event):
            return len(self._event_source_paths(event)) == exact_support

        other_filtered = [event for event in other_events if has_exact_support(event)]
        tra_filtered = [event for event in tra_events if has_exact_support(event)]
        bnd_filtered = [event for event in bnd_events if has_exact_support(event)]
        return other_filtered + tra_filtered + bnd_filtered

    def get_events_by_support_range(self, min_support=None, max_support=None):
        """Get events supported by a range of distinct physical input sources."""
        other_events, tra_events, bnd_events = self._all_selectable_events()

        def within_range(event):
            support_count = len(self._event_source_paths(event))
            if min_support is not None and support_count < min_support:
                return False
            return not (max_support is not None and support_count > max_support)

        other_filtered = [event for event in other_events if within_range(event)]
        tra_filtered = [event for event in tra_events if within_range(event)]
        bnd_filtered = [event for event in bnd_events if within_range(event)]
        return other_filtered + tra_filtered + bnd_filtered

    @staticmethod
    def _expression_identifier(file_path):
        return _expression_identifier(file_path)

    def _expression_source_map(self):
        """Build unambiguous basename/identifier mappings for --expression."""
        return _build_expression_source_entries(
            self.all_input_files,
            normalize_func=self._normalize_source,
        )

    def _prepare_expression(self, expression):
        """Prepare one expression once for repeated event evaluation."""
        entries = self._expression_source_map()
        compiled, referenced_identifiers = _compile_expression(expression, entries)
        referenced_entries = tuple(
            entry
            for entry in entries
            if entry[1] in referenced_identifiers
        )
        return compiled, referenced_entries

    @staticmethod
    def _evaluate_prepared_expression(compiled, entries, normalized_event_sources):
        context = {
            identifier: normalized in normalized_event_sources
            for _basename, identifier, normalized in entries
        }
        try:
            return bool(eval(compiled, {"__builtins__": {}}, context))
        except Exception as exc:
            raise ValueError(f"Invalid expression: {exc}") from exc

    def get_events_by_expression(self, expression):
        """Get events that satisfy a logical source expression."""
        compiled, entries = self._prepare_expression(expression)
        other_events, tra_events, bnd_events = self._all_selectable_events()

        def selected(event):
            return self._evaluate_prepared_expression(
                compiled,
                entries,
                self._event_source_paths(event),
            )

        other_filtered = [event for event in other_events if selected(event)]
        tra_filtered = [event for event in tra_events if selected(event)]
        bnd_filtered = [event for event in bnd_events if selected(event)]
        return other_filtered + tra_filtered + bnd_filtered

    def evaluate_expression(self, expression, event_sources):
        """Evaluate one logical expression against event-source identities."""
        compiled, entries = self._prepare_expression(expression)
        normalized_event_sources = {
            self._normalize_source(source)
            for source in event_sources
        }
        return self._evaluate_prepared_expression(
            compiled,
            entries,
            normalized_event_sources,
        )
