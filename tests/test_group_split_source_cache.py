"""Source-path caching must not change the #192 group-splitting decision."""

from types import SimpleNamespace

import octopusv.merger.merge_split as merge_split
from octopusv.merger.merge_split import GroupSplitter


def _del(source, length, *, original=None):
    event = SimpleNamespace(
        source_file=str(source),
        sv_type="DEL",
        info={"SVTYPE": "DEL", "SVLEN": str(length)},
        start_pos=1000,
        end_pos=1000 + length,
    )
    if original is not None:
        event._original_source_file = str(original)
    return event


def test_source_path_is_normalized_once_per_distinct_spelling(tmp_path, monkeypatch):
    caller_a = tmp_path / "caller_a.svcf"
    caller_b = tmp_path / "caller_b.svcf"
    alias_a = tmp_path / "alias_a.svcf"
    caller_a.touch()
    caller_b.touch()
    alias_a.symlink_to(caller_a)

    original_normalize = merge_split.normalize_source_path
    normalized = []

    def counted(path):
        normalized.append(str(path))
        return original_normalize(path)

    monkeypatch.setattr(merge_split, "normalize_source_path", counted)
    splitter = GroupSplitter("caller")

    for _ in range(20):
        group = [_del(caller_a, 100), _del(alias_a, 200), _del(caller_b, 101)]
        parts = splitter.split(group)
        # Aliased files are the SAME caller. Incompatible lengths must split.
        assert [[member.info["SVLEN"] for member in part] for part in parts] == [
            ["100", "101"], ["200"],
        ]

    assert normalized == [str(caller_a), str(alias_a), str(caller_b)]


def test_source_cache_lifetime_is_one_splitter_instance(tmp_path):
    original = tmp_path / "original.svcf"
    changed = tmp_path / "changed.svcf"
    alias = tmp_path / "alias.svcf"
    original.touch()
    changed.touch()
    alias.symlink_to(original)

    initial = GroupSplitter("caller")
    assert initial._source_key(_del(alias, 100)) == initial._source_key(_del(original, 100))

    alias.unlink()
    alias.symlink_to(changed)
    new_run = GroupSplitter("caller")
    assert new_run._source_key(_del(alias, 100)) == new_run._source_key(_del(changed, 100))
    assert new_run._source_key(_del(alias, 100)) != new_run._source_key(_del(original, 100))


def test_original_source_path_overrides_mutated_representative(tmp_path):
    a = tmp_path / "a.svcf"
    b = tmp_path / "b.svcf"
    a.touch()
    b.touch()
    splitter = GroupSplitter("caller")
    event = _del("a.svcf,b.svcf", 100, original=a)
    assert splitter._source_key(event) == merge_split.normalize_source_path(a)


def test_sample_mode_still_checks_direct_match(tmp_path):
    a = tmp_path / "a.svcf"
    a.touch()
    same_source = [_del(a, 100), _del(a, 102)]
    assert len(GroupSplitter("caller").split(same_source)) == 1
    sample_splitter = GroupSplitter("sample", direct_match=lambda x, y: False)
    assert len(sample_splitter.split(same_source)) == 2
