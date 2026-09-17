from types import SimpleNamespace

from octopusv.merger.name_mapper import NameMapper
from octopusv.merger.sv_merge_writer import MergeWriterMixin


FORMAT = "GT:AD:LN:ST:QV:TY:ID:SC:REF:ALT:CO"


class DummyWriter(MergeWriterMixin):
    def __init__(self, input_files):
        self.all_input_files = [str(path) for path in input_files]


def _block(gt, record_id, caller):
    return (
        f"{gt}:5,5:50:+-:60:INS:{record_id}:{caller}:N:<INS>:"
        "chr1_100-chr1_150"
    )


def _sample_payload(record_id, sources, blocks):
    return {
        "ID": blocks[0].split(":")[6],
        "GT": blocks[0].split(":")[0],
        "_octopusv_evidence_payload": {
            "format": FORMAT,
            "blocks": tuple(blocks),
            "sources": ",".join(sources),
            "source_ids": ",".join(
                block.split(":")[6] for block in blocks
            ),
            "record": {
                "id": record_id,
                "ref": "N",
                "alt": "<INS>",
                "quality": "60",
                "svtype": "INS",
                "strand": "+-",
                "svlen": "50",
                "chrom": "chr1",
                "pos": 100,
                "end_chrom": "chr1",
                "end_pos": 150,
            },
        },
    }


def _event(records):
    merged_samples = [
        (sample_name, sample_format, sample_data)
        for _, sample_name, sample_format, sample_data in records
    ]
    return SimpleNamespace(
        sv_id="merged.1",
        source_file=",".join(str(record[0]) for record in records),
        merged_samples=merged_samples,
        merged_sample_records=records,
    )


def test_sample_mode_synthesizes_order_independent_calls(tmp_path):
    sample_a = tmp_path / "sampleA.svcf"
    sample_b = tmp_path / "sampleB.svcf"
    writer = DummyWriter([sample_a, sample_b])
    mapper = NameMapper(
        [str(sample_a), str(sample_b)],
        mode="sample",
        custom_names=["A", "B"],
    )

    data_a = _sample_payload(
        "A.record",
        ["cuteSV", "svim", "pbsv"],
        [
            _block("0/1", "cute.a", "cuteSV"),
            _block("0/1", "svim.a", "svim"),
            _block("1/1", "pbsv.a", "pbsv"),
        ],
    )
    data_b = _sample_payload(
        "B.record",
        ["pbsv", "svim", "cuteSV"],
        [
            _block("1/1", "pbsv.b", "pbsv"),
            _block("0/1", "svim.b", "svim"),
            _block("0/1", "cute.b", "cuteSV"),
        ],
    )

    event = _event(
        [
            (str(sample_a), "SAMPLE", FORMAT, data_a),
            (str(sample_b), "SAMPLE", FORMAT, data_b),
        ]
    )

    processed, summary = writer._prepare_events_for_sample_mode([event], mapper)
    sample_a_out, sample_b_out = processed[0].ordered_samples

    assert sample_a_out["GT"] == "1/."
    assert sample_b_out["GT"] == "1/."
    assert sample_a_out["UC"] == sample_b_out["UC"] == "3"
    assert sample_a_out["UV"] == sample_b_out["UV"] == "3"
    assert sample_a_out["AD"] == sample_b_out["AD"] == ".,."
    assert sample_a_out["SC"] == sample_b_out["SC"] == "OctopuSV"
    assert sample_a_out["ID"] == "A.record"
    assert sample_b_out["ID"] == "B.record"
    assert summary["legacy_sample_events"] == 0


def test_sample_mode_duplicate_source_contributes_one_caller_state(tmp_path):
    sample_a = tmp_path / "sampleA.svcf"
    writer = DummyWriter([sample_a])
    mapper = NameMapper([str(sample_a)], mode="sample", custom_names=["A"])

    data = _sample_payload(
        "A.dup",
        ["cuteSV", "cuteSV", "svim"],
        [
            _block("0/1", "cute.1", "cuteSV"),
            _block("1/1", "cute.2", "cuteSV"),
            _block("0/1", "svim.1", "svim"),
        ],
    )
    event = _event([(str(sample_a), "SAMPLE", FORMAT, data)])

    processed, _ = writer._prepare_events_for_sample_mode([event], mapper)
    sample = processed[0].ordered_samples[0]

    assert sample["GT"] == "1/."
    assert sample["UC"] == "2"
    assert sample["UV"] == "2"


def test_sample_mode_combines_all_records_from_same_input_for_consensus(tmp_path):
    sample_a = tmp_path / "sampleA.svcf"
    writer = DummyWriter([sample_a])
    mapper = NameMapper([str(sample_a)], mode="sample", custom_names=["A"])

    first = _sample_payload(
        "A.first",
        ["cuteSV"],
        [_block("0/1", "cute.1", "cuteSV")],
    )
    second = _sample_payload(
        "A.second",
        ["svim"],
        [_block("1/1", "svim.1", "svim")],
    )
    event = _event(
        [
            (str(sample_a), "SAMPLE", FORMAT, first),
            (str(sample_a), "SAMPLE", FORMAT, second),
        ]
    )

    processed, _ = writer._prepare_events_for_sample_mode([event], mapper)
    sample = processed[0].ordered_samples[0]

    assert sample["GT"] == "1/."
    assert sample["UC"] == "2"
    assert sample["UV"] == "2"
    # Representation fields come from the first deterministic input record;
    # genotype consensus uses evidence from both records.
    assert sample["ID"] == "A.first"


def test_legacy_caller_unresolved_evidence_is_counted_not_silent(tmp_path):
    input_file = tmp_path / "callerA.svcf"
    writer = DummyWriter([input_file])

    event = SimpleNamespace(
        sv_id="legacy.1",
        source_file="does_not_match_any_input",
        merged_samples=[
            ("SAMPLE", "FMT", {"ID": "legacy.id"}),
        ],
    )

    records = writer._prepare_caller_records_legacy(event)

    assert records == []
    assert writer._legacy_caller_event_count == 1
    assert writer._legacy_caller_unresolved_evidence_count == 1


def test_sample_mode_synthesized_fields_use_record_level_payload(tmp_path):
    sample_a = tmp_path / "sampleA.svcf"
    writer = DummyWriter([sample_a])
    mapper = NameMapper([str(sample_a)], mode="sample", custom_names=["A"])

    data = _sample_payload(
        "A.record.level.id",
        ["cuteSV"],
        [_block("0/1", "caller.level.id", "cuteSV")],
    )
    payload = data["_octopusv_evidence_payload"]
    payload["record"].update(
        {
            "ref": "G",
            "alt": "<DEL>",
            "quality": "42",
            "svtype": "DEL",
            "strand": "-+",
            "svlen": "75",
            "chrom": "chr2",
            "pos": 200,
            "end_chrom": "chr2",
            "end_pos": 275,
        }
    )
    event = _event([(str(sample_a), "SAMPLE", FORMAT, data)])

    processed, _ = writer._prepare_events_for_sample_mode([event], mapper)
    sample = processed[0].ordered_samples[0]

    assert sample["ID"] == "A.record.level.id"
    assert sample["REF"] == "G"
    assert sample["ALT"] == "<DEL>"
    assert sample["QV"] == "42"
    assert sample["TY"] == "DEL"
    assert sample["ST"] == "-+"
    assert sample["LN"] == "75"
    assert sample["CO"] == "chr2_200-chr2_275"


def test_sample_mode_refuses_to_guess_missing_caller_identity(tmp_path):
    sample_a = tmp_path / "sampleA.svcf"
    writer = DummyWriter([sample_a])
    mapper = NameMapper([str(sample_a)], mode="sample", custom_names=["A"])

    block = _block("0/1", "unknown.1", ".")
    data = _sample_payload("A.record", ["."], [block])
    event = _event([(str(sample_a), "SAMPLE", FORMAT, data)])

    import pytest

    with pytest.raises(ValueError, match="no explicit caller identity"):
        writer._prepare_events_for_sample_mode([event], mapper)
