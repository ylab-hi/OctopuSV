from __future__ import annotations

from types import SimpleNamespace

from octopusv.merger.multi_sample_writer import MultiSampleWriter
from octopusv.merger.name_mapper import NameMapper
from octopusv.merger.sv_merge_writer import MergeWriterMixin
from octopusv.utils.vcf_info import format_vcf_info_item


CALLER_FORMAT = "GT:AD:LN:ST:QV:TY:ID:SC:REF:ALT:CO"


def _caller_sample_data():
    return {
        "GT": "0/1",
        "AD": "5,5",
        "LN": "100",
        "ST": ".",
        "QV": "60",
        "TY": "DEL",
        "ID": "caller.record.1",
        "SC": "cuteSV",
        "REF": "N",
        "ALT": "<DEL>",
        "CO": "chr1_100-chr1_200",
    }


def _sample_mode_data():
    return {
        "GT": "0/1",
        "AD": ".,.",
        "UC": "1",
        "UV": "1",
        "LN": "100",
        "ST": ".",
        "QV": "60",
        "TY": "DEL",
        "ID": "sample.event.1",
        "SC": "OctopuSV",
        "REF": "N",
        "ALT": "<DEL>",
        "CO": "chr1_100-chr1_200",
    }


def _event(*, sample_data, source_file="/inputs/cute.svcf"):
    return SimpleNamespace(
        chrom="chr1",
        pos=100,
        sv_id="merged.1",
        ref="N",
        alt="<DEL>",
        quality="60",
        filter="PASS",
        info={
            "SVTYPE": "DEL",
            "END": "200",
            "SVLEN": "-100",
            "PRECISE": True,
            "LITERAL_TRUE": "True",
        },
        format=CALLER_FORMAT,
        sample=sample_data,
        source_file=source_file,
        merged_samples=[("SAMPLE", CALLER_FORMAT, sample_data)],
        merged_sample_records=[
            (source_file, "SAMPLE", CALLER_FORMAT, sample_data)
        ],
    )


class _CallerWriter(MergeWriterMixin):
    def __init__(self, input_files):
        self.all_input_files = [str(path) for path in input_files]


def _record_info_field(path):
    record_lines = [
        line
        for line in path.read_text().splitlines()
        if line and not line.startswith("#")
    ]
    assert len(record_lines) == 1
    return record_lines[0].split("\t")[7]


def test_format_vcf_info_item_only_treats_boolean_true_as_flag():
    assert format_vcf_info_item("PRECISE", True) == "PRECISE"
    assert format_vcf_info_item("LITERAL_TRUE", "True") == "LITERAL_TRUE=True"
    assert format_vcf_info_item("COUNT", 3) == "COUNT=3"


def test_caller_writer_serializes_info_flag_as_bare_key(monkeypatch, tmp_path):
    input_file = tmp_path / "caller.svcf"
    input_file.write_text("")

    writer = _CallerWriter([input_file])
    monkeypatch.setattr(
        writer,
        "_write_vcf_header",
        lambda handle, contigs, input_files: None,
    )

    output = tmp_path / "caller_merged.svcf"
    writer.write_results(
        output,
        [_event(sample_data=_caller_sample_data(), source_file=str(input_file))],
        {},
        mode="caller",
        input_files=[str(input_file)],
    )

    info = _record_info_field(output)
    tokens = info.split(";")

    assert "PRECISE" in tokens
    assert "PRECISE=True" not in tokens
    assert "LITERAL_TRUE=True" in tokens


def test_sample_writer_serializes_info_flag_as_bare_key(tmp_path):
    input_file = tmp_path / "sample.svcf"
    mapper = NameMapper([str(input_file)], mode="sample")
    writer = MultiSampleWriter(mapper)

    event = _event(sample_data=_sample_mode_data(), source_file=str(input_file))
    event.ordered_samples = [_sample_mode_data()]

    output = tmp_path / "sample_merged.svcf"
    writer.write_results(output, [event], {"chr1": 1000000})

    info = _record_info_field(output)
    tokens = info.split(";")

    assert "PRECISE" in tokens
    assert "PRECISE=True" not in tokens
    assert "LITERAL_TRUE=True" in tokens
