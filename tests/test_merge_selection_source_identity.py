from types import SimpleNamespace

import pytest

from octopusv.merger.sv_merge_selection import MergeSelectionMixin


class _EmptySpecialMerger:
    def get_merged_events(self):
        return []


class _SelectionHarness(MergeSelectionMixin):
    def __init__(self, input_files, events):
        self.all_input_files = [str(path) for path in input_files]
        self._events = list(events)
        self.tra_merger = _EmptySpecialMerger()
        self.bnd_merger = _EmptySpecialMerger()

    def get_all_merged_events(self):
        return list(self._events)


def _event(name, *sources):
    return SimpleNamespace(
        sv_id=name,
        sv_type="DEL",
        source_file=",".join(str(source) for source in sources),
    )


def _ids(events):
    return {event.sv_id for event in events}


def test_union_distinguishes_same_basename_in_different_directories(tmp_path):
    source_a = tmp_path / "d1" / "shared.svcf"
    source_b = tmp_path / "d2" / "shared.svcf"

    merger = _SelectionHarness(
        [source_a, source_b],
        [
            _event("A", source_a),
            _event("B", source_b),
            _event("AB", source_a, source_b),
        ],
    )

    selected = merger.get_events_by_source([source_a], operation="union")

    assert _ids(selected) == {"A", "AB"}


def test_intersection_distinguishes_same_basename_in_different_directories(tmp_path):
    source_a = tmp_path / "d1" / "shared.svcf"
    source_b = tmp_path / "d2" / "shared.svcf"

    merger = _SelectionHarness(
        [source_a, source_b],
        [
            _event("A", source_a),
            _event("B", source_b),
            _event("AB", source_a, source_b),
        ],
    )

    selected = merger.get_events_by_source(
        [source_a, source_b],
        operation="intersection",
    )

    assert _ids(selected) == {"AB"}


def test_specific_distinguishes_same_basename_in_different_directories(tmp_path):
    source_a = tmp_path / "d1" / "shared.svcf"
    source_b = tmp_path / "d2" / "shared.svcf"

    merger = _SelectionHarness(
        [source_a, source_b],
        [
            _event("A", source_a),
            _event("B", source_b),
            _event("AB", source_a, source_b),
        ],
    )

    selected = merger.get_events_by_source([source_a], operation="specific")

    assert _ids(selected) == {"A"}


def test_specific_multiple_sources_keeps_events_supported_only_by_requested_set(tmp_path):
    source_a = tmp_path / "a.svcf"
    source_b = tmp_path / "b.svcf"
    source_c = tmp_path / "c.svcf"

    merger = _SelectionHarness(
        [source_a, source_b, source_c],
        [
            _event("A", source_a),
            _event("B", source_b),
            _event("C", source_c),
            _event("AB", source_a, source_b),
            _event("AC", source_a, source_c),
            _event("ABC", source_a, source_b, source_c),
        ],
    )

    selected = merger.get_events_by_source(
        [source_a, source_b],
        operation="specific",
    )

    assert _ids(selected) == {"A", "B", "AB"}


def test_specific_rejects_source_that_is_not_a_merge_input(tmp_path):
    source_a = tmp_path / "a.svcf"
    unknown = tmp_path / "unknown.svcf"
    merger = _SelectionHarness([source_a], [_event("A", source_a)])

    with pytest.raises(ValueError, match="not one of the merge inputs"):
        merger.get_events_by_source([unknown], operation="specific")


def test_source_matching_uses_realpath_identity(tmp_path):
    real_source = tmp_path / "real.svcf"
    real_source.write_text("")
    alias = tmp_path / "alias.svcf"
    alias.symlink_to(real_source)

    merger = _SelectionHarness(
        [real_source],
        [_event("A", alias)],
    )

    selected = merger.get_events_by_source([real_source], operation="union")

    assert _ids(selected) == {"A"}


def test_expression_rejects_same_basename_from_different_directories(tmp_path):
    source_a = tmp_path / "d1" / "shared.svcf"
    source_b = tmp_path / "d2" / "shared.svcf"
    merger = _SelectionHarness(
        [source_a, source_b],
        [_event("A", source_a), _event("B", source_b)],
    )

    with pytest.raises(ValueError, match="share the basename"):
        merger.get_events_by_expression("shared.svcf")


def test_expression_rejects_sanitized_identifier_collision(tmp_path):
    source_a = tmp_path / "a-b.svcf"
    source_b = tmp_path / "a_b.svcf"
    merger = _SelectionHarness(
        [source_a, source_b],
        [_event("A", source_a), _event("B", source_b)],
    )

    with pytest.raises(ValueError, match="same expression identifier"):
        merger.get_events_by_expression("a-b.svcf OR a_b.svcf")


def test_expression_uses_exact_source_identity_for_unique_basenames(tmp_path):
    source_a = tmp_path / "a.svcf"
    source_b = tmp_path / "b.svcf"

    merger = _SelectionHarness(
        [source_a, source_b],
        [
            _event("A", source_a),
            _event("B", source_b),
            _event("AB", source_a, source_b),
        ],
    )

    selected = merger.get_events_by_expression("a.svcf AND NOT b.svcf")

    assert _ids(selected) == {"A"}


def test_exact_support_uses_normalized_source_identity(tmp_path):
    real_source = tmp_path / "real.svcf"
    real_source.write_text("")
    alias = tmp_path / "alias.svcf"
    alias.symlink_to(real_source)

    merger = _SelectionHarness(
        [real_source],
        [_event("A", real_source, alias)],
    )

    selected = merger.get_events_by_exact_support(1)

    assert _ids(selected) == {"A"}


def test_expression_normalization_is_cached_across_many_events(tmp_path, monkeypatch):
    import octopusv.merger.sv_merge_selection as selection_module

    sources = [tmp_path / f"source_{idx}.svcf" for idx in range(20)]
    events = [
        _event(f"event_{idx}", sources[idx % len(sources)])
        for idx in range(500)
    ]
    merger = _SelectionHarness(sources, events)

    original = selection_module.normalize_source_path
    calls = []

    def counted(source):
        calls.append(str(source))
        return original(source)

    monkeypatch.setattr(selection_module, "normalize_source_path", counted)

    selected = merger.get_events_by_expression(
        "source_0.svcf OR source_1.svcf"
    )

    assert len(selected) == 50
    # Each raw input source string is normalized at most once. The old
    # implementation normalized all inputs again for every event.
    assert len(calls) == len(sources)
    assert set(calls) == {str(source) for source in sources}


def test_expression_unknown_identifier_fails_before_event_evaluation(tmp_path):
    source_a = tmp_path / "a.svcf"
    merger = _SelectionHarness([source_a], [_event("A", source_a)])

    with pytest.raises(ValueError, match="unknown source identifier"):
        merger.get_events_by_expression("a.svcf AND missing")


@pytest.mark.parametrize(
    "expression",
    [
        "a.svcf + b.svcf",
        "a.svcf == b.svcf",
        "a.svcf.__class__",
    ],
)
def test_expression_rejects_non_boolean_ast_nodes(tmp_path, expression):
    source_a = tmp_path / "a.svcf"
    source_b = tmp_path / "b.svcf"
    merger = _SelectionHarness(
        [source_a, source_b],
        [_event("A", source_a), _event("B", source_b)],
    )

    with pytest.raises(
        ValueError,
        match="only input source names combined with AND, OR, NOT",
    ):
        merger.get_events_by_expression(expression)
