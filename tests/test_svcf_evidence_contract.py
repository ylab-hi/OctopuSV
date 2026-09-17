from octopusv.merger.multi_sample_writer import MultiSampleWriter
from octopusv.merger.name_mapper import NameMapper
from octopusv.merger.sv_merge_writer import MergeWriterMixin


FORMAT_FIELD = "GT:AD:LN:ST:QV:TY:ID:SC:REF:ALT:CO"
FORMAT_KEYS = FORMAT_FIELD.split(":")


def _sample_data(**overrides):
    data = {
        "GT": "0/1",
        "AD": "5,5",
        "LN": "50",
        "ST": ".",
        "QV": "60",
        "TY": "INS",
        "ID": "evidence.1",
        "SC": "cuteSV",
        "REF": "N",
        "ALT": ".",
        "CO": ".",
    }
    data.update(overrides)
    return data


def test_caller_writer_keeps_full_fixed_width_evidence_block():
    text = MergeWriterMixin().format_sample_values(
        FORMAT_KEYS,
        _sample_data(),
    )

    assert len(text.split(":")) == len(FORMAT_KEYS)
    assert text.endswith(":N:.:.")


def test_source_id_prefers_evidence_id_over_record_level_fallback():
    sample_data = {
        "ID": "evidence.id",
        "original_id": "representative.record.id",
    }

    assert (
        MergeWriterMixin._source_id_from_sample_data(sample_data)
        == "evidence.id"
    )


def test_source_id_uses_original_id_only_when_evidence_id_absent():
    assert (
        MergeWriterMixin._source_id_from_sample_data(
            {"original_id": "legacy.record.id"}
        )
        == "legacy.record.id"
    )


def test_explicit_missing_evidence_id_stays_missing():
    sample_data = {
        "ID": ".",
        "original_id": "representative.record.id",
    }

    assert MergeWriterMixin._source_id_from_sample_data(sample_data) == "."


def test_sample_writer_prefers_evidence_id():
    mapper = NameMapper(
        ["/tmp/sample1.svcf"],
        mode="sample",
    )
    writer = MultiSampleWriter(mapper)

    sample_data = {
        "ID": "evidence.id",
        "original_id": "representative.record.id",
    }

    assert (
        writer._sample_id_from_data(sample_data, FORMAT_KEYS)
        == "evidence.id"
    )


def test_sample_writer_keeps_all_format_fields_for_real_and_missing_columns():
    mapper = NameMapper(
        ["/tmp/sample1.svcf", "/tmp/sample2.svcf"],
        mode="sample",
    )
    writer = MultiSampleWriter(mapper)

    columns = writer._format_sample_columns(
        [_sample_data(), None],
        FORMAT_KEYS,
    ).split("\t")

    assert len(columns) == 2
    assert len(columns[0].split(":")) == len(FORMAT_KEYS)
    assert len(columns[1].split(":")) == len(FORMAT_KEYS)

def test_sample_writer_normalizes_none_and_list_values():
    mapper = NameMapper(
        ["/tmp/sample1.svcf"],
        mode="sample",
    )
    writer = MultiSampleWriter(mapper)

    sample_data = _sample_data(AD=[5, 7], QV=None)
    column = writer._format_sample_columns(
        [sample_data],
        FORMAT_KEYS,
    )

    fields = column.split(":")
    assert len(fields) == len(FORMAT_KEYS)
    assert fields[1] == "5,7"
    assert fields[4] == "."

