"""Regression guards for caller-mode write-scoped source path caching."""
from types import SimpleNamespace

import pytest

from octopusv.merger.sv_merge_writer import MergeWriterMixin


class Writer(MergeWriterMixin):
    pass


def _event(pos, records):
    evidence = [
        (path, "SAMPLE", "GT:ID:SC", {"GT": "0/1", "ID": record_id, "SC": caller})
        for path, record_id, caller in records
    ]
    return SimpleNamespace(
        chrom="chr1", pos=pos, sv_id=f"event.{pos}", ref="N", alt="<DEL>",
        quality="60", filter="PASS", info={"SVTYPE": "DEL", "END": str(pos + 50)},
        format="GT:ID:SC", merged_sample_records=evidence, source_file="",
    )


def _data_rows(path):
    return [line.split("\t") for line in path.read_text().splitlines() if not line.startswith("#")]


def test_caller_cache_preserves_symlink_mapping_and_evidence_order(monkeypatch, tmp_path):
    a = tmp_path / "alpha.svcf"
    b = tmp_path / "beta.svcf"
    alias = tmp_path / "alpha_alias.svcf"
    a.touch()
    b.touch()
    alias.symlink_to(a)

    writer = Writer()
    writer.all_input_files = [str(a), str(b)]
    monkeypatch.setattr(writer, "_write_vcf_header", lambda *_: None)

    original_normalize = writer._normalize_source_path
    calls = []

    def counted(path):
        calls.append(str(path))
        return original_normalize(path)

    monkeypatch.setattr(writer, "_normalize_source_path", counted)
    events = [
        _event(pos, [(str(b), f"B.{pos}", "beta"), (str(alias), f"A.{pos}", "alpha")])
        for pos in range(100, 1100)
    ]
    out = tmp_path / "merged.svcf"
    writer.write_results(out, events, {"chr1": 1000000}, mode="caller", input_files=[str(a), str(b)])

    # Two input paths + two distinct evidence-source paths, only once per write.
    assert len(calls) == 4
    rows = _data_rows(out)
    assert len(rows) == 1000
    for row in rows:
        pos = int(row[1])
        assert "SOURCES=alpha,beta" in row[7]
        assert f"SOURCE_IDS=A.{pos},B.{pos}" in row[7]
        assert [block.split(":")[1] for block in row[9:]] == [f"A.{pos}", f"B.{pos}"]


def test_caller_cache_keeps_duplicate_realpath_ambiguity(tmp_path, monkeypatch):
    a = tmp_path / "alpha.svcf"
    alias = tmp_path / "alias.svcf"
    a.touch()
    alias.symlink_to(a)

    writer = Writer()
    writer.all_input_files = [str(a), str(alias)]
    monkeypatch.setattr(writer, "_write_vcf_header", lambda *_: None)

    event = _event(100, [(str(a), "evidence.1", "alpha")])
    with pytest.raises(ValueError, match="matches more than one merge input"):
        writer.write_results(
            tmp_path / "out.svcf", [event], {"chr1": 1000000},
            mode="caller", input_files=[str(a), str(alias)],
        )
