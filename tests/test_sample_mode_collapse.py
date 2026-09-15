from types import SimpleNamespace

from octopusv.cli.merge import _mark_sample_input_collapses
from octopusv.merger.name_mapper import NameMapper
from octopusv.merger.sv_merge_writer import MergeWriterMixin


class DummyWriter(MergeWriterMixin):
    def __init__(self, input_files):
        self.all_input_files = [str(path) for path in input_files]


def _sample(record_id, collapsed=0):
    data = {
        "GT": "0/1",
        "AD": "5,5",
        "LN": "50",
        "ST": ".",
        "QV": "60",
        "TY": "INS",
        "ID": record_id,
        "SC": "caller",
        "REF": "N",
        "ALT": "<INS>",
        "CO": "chr1_100-chr1_150",
    }
    if collapsed:
        data["_octopusv_collapsed_evidence_count"] = collapsed
    return data


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


def test_sample_mode_keeps_first_deterministic_evidence_and_counts_extra(tmp_path):
    sample_a = tmp_path / "sampleA.svcf"
    sample_b = tmp_path / "sampleB.svcf"
    writer = DummyWriter([sample_a, sample_b])
    mapper = NameMapper(
        [str(sample_a), str(sample_b)],
        mode="sample",
        custom_names=["A", "B"],
    )

    first_a = _sample("a.first")
    second_a = _sample("a.second")
    first_b = _sample("b.first")
    event = _event(
        [
            (str(sample_a), "SAMPLE", "FMT", first_a),
            (str(sample_a), "SAMPLE", "FMT", second_a),
            (str(sample_b), "SAMPLE", "FMT", first_b),
        ]
    )

    processed, summary = writer._prepare_events_for_sample_mode([event], mapper)

    assert processed[0].ordered_samples == [first_a, first_b]
    assert summary["sample_collapse_records"] == 1
    assert summary["sample_collapsed_evidence_blocks"] == 1


def test_sample_mode_counts_nested_caller_evidence_on_retained_sample(tmp_path):
    sample_a = tmp_path / "sampleA.svcf"
    writer = DummyWriter([sample_a])
    mapper = NameMapper(
        [str(sample_a)],
        mode="sample",
        custom_names=["A"],
    )

    collapsed_input = _sample("callerA.first", collapsed=2)
    event = _event(
        [(str(sample_a), "SAMPLE", "FMT", collapsed_input)]
    )

    processed, summary = writer._prepare_events_for_sample_mode([event], mapper)

    assert processed[0].ordered_samples == [collapsed_input]
    assert summary["sample_collapse_records"] == 1
    assert summary["sample_collapsed_evidence_blocks"] == 2


def test_sample_mode_counts_all_underlying_blocks_when_duplicate_record_is_collapsed(tmp_path):
    sample_a = tmp_path / "sampleA.svcf"
    writer = DummyWriter([sample_a])
    mapper = NameMapper(
        [str(sample_a)],
        mode="sample",
        custom_names=["A"],
    )

    retained = _sample("a.first", collapsed=1)
    dropped = _sample("a.second", collapsed=2)
    event = _event(
        [
            (str(sample_a), "SAMPLE", "FMT", retained),
            (str(sample_a), "SAMPLE", "FMT", dropped),
        ]
    )

    _, summary = writer._prepare_events_for_sample_mode([event], mapper)

    # retained block represents one additional nested block; the second record
    # contributes its selected block plus two nested blocks = three more.
    assert summary["sample_collapse_records"] == 1
    assert summary["sample_collapsed_evidence_blocks"] == 4


def test_legacy_caller_unresolved_evidence_is_counted_not_silent(tmp_path):
    input_file = tmp_path / "callerA.svcf"
    writer = DummyWriter([input_file])

    event = SimpleNamespace(
        sv_id="legacy.1",
        source_file="does_not_match_any_input",
        merged_samples=[
            ("SAMPLE", "FMT", _sample("legacy.id")),
        ],
    )

    records = writer._prepare_caller_records_legacy(event)

    assert records == []
    assert writer._legacy_caller_event_count == 1
    assert writer._legacy_caller_unresolved_evidence_count == 1


def test_sample_mode_retains_each_inputs_first_block_when_caller_order_differs(tmp_path):
    """Different per-sample caller orders remain deterministic and visible."""
    sample_a = tmp_path / "sampleA.svcf"
    sample_b = tmp_path / "sampleB.svcf"
    writer = DummyWriter([sample_a, sample_b])
    mapper = NameMapper(
        [str(sample_a), str(sample_b)],
        mode="sample",
        custom_names=["A", "B"],
    )

    # sample A was caller-merged as cuteSV -> SVIM; sample B used the reverse
    # order. SVCFEvent.sample represents the first evidence column from each
    # input, and _mark_sample_input_collapses records the additional block.
    first_a = _sample("cuteSV.a")
    first_b = _sample("SVIM.b")
    input_event_a = SimpleNamespace(
        sample=first_a,
        raw_sample_columns=["cuteSV-block", "SVIM-block"],
    )
    input_event_b = SimpleNamespace(
        sample=first_b,
        raw_sample_columns=["SVIM-block", "cuteSV-block"],
    )
    _mark_sample_input_collapses([input_event_a, input_event_b])

    merged_event = _event(
        [
            (str(sample_a), "SAMPLE", "FMT", first_a),
            (str(sample_b), "SAMPLE", "FMT", first_b),
        ]
    )

    processed, summary = writer._prepare_events_for_sample_mode(
        [merged_event],
        mapper,
    )

    assert [sample["ID"] for sample in processed[0].ordered_samples] == [
        "cuteSV.a",
        "SVIM.b",
    ]
    assert summary["sample_collapse_records"] == 1
    assert summary["sample_collapsed_evidence_blocks"] == 2
