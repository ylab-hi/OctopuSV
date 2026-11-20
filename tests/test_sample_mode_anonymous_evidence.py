"""Sample mode with evidence that names no caller (no SOURCES, SC=".").

A caller that writes no ##source (e.g. SVDSS) gives `correct` output with
SC="." on every record.  Such a single-caller input must merge in sample mode:
its records are one voter, exactly as if SC named the caller.  Evidence with
and without a caller name in one sample of one merged event is ambiguous and
is refused rather than counted as separate callers.
"""

from types import SimpleNamespace

import pytest
from typer.testing import CliRunner

from octopusv.cli.cli import app
from octopusv.merger.name_mapper import NameMapper
from octopusv.merger.sv_merge_writer import MergeWriterMixin

FORMAT = "GT:AD:LN:ST:QV:TY:ID:SC:REF:ALT:CO"


class DummyWriter(MergeWriterMixin):
    def __init__(self, input_files):
        self.all_input_files = [str(path) for path in input_files]


def _block(gt, record_id, caller):
    return f"{gt}:5,5:50:+-:60:INS:{record_id}:{caller}:N:<INS>:chr1_100-chr1_150"


def _sample_payload(record_id, sources, blocks):
    return {
        "ID": blocks[0].split(":")[6],
        "GT": blocks[0].split(":")[0],
        "_octopusv_evidence_payload": {
            "format": FORMAT,
            "blocks": tuple(blocks),
            "sources": ",".join(sources),
            "source_ids": ",".join(block.split(":")[6] for block in blocks),
            "record": {
                "id": record_id, "ref": "N", "alt": "<INS>", "quality": "60",
                "svtype": "INS", "strand": "+-", "svlen": "50",
                "chrom": "chr1", "pos": 100, "end_chrom": "chr1", "end_pos": 150,
            },
        },
    }


def _event(records):
    return SimpleNamespace(
        sv_id="merged.1",
        source_file=",".join(str(record[0]) for record in records),
        merged_samples=[record[1:] for record in records],
        merged_sample_records=records,
    )


def _prepare(tmp_path, records):
    """records: (record_id, GT, SC) single-block records of one input file."""
    sample_a = tmp_path / "sampleA.svcf"
    writer = DummyWriter([sample_a])
    mapper = NameMapper([str(sample_a)], mode="sample", custom_names=["A"])
    event = _event([
        (
            str(sample_a),
            "SAMPLE",
            FORMAT,
            _sample_payload(record_id, ["."], [_block(gt, record_id, caller)]),
        )
        for record_id, gt, caller in records
    ])
    processed, _ = writer._prepare_events_for_sample_mode([event], mapper)
    return processed[0].ordered_samples[0]


def test_unnamed_single_caller_records_are_one_voter(tmp_path):
    unnamed = _prepare(tmp_path, [("A.r1", "0/1", "."), ("A.r2", "1/1", ".")])
    named = _prepare(tmp_path, [("A.r1", "0/1", "SVDSS"), ("A.r2", "1/1", "SVDSS")])

    assert unnamed["UV"] == "1"
    assert unnamed["UC"] == "1"
    assert unnamed == named


def test_lone_unnamed_record_matches_named_record(tmp_path):
    for gt in ("0/1", "1/1", "0/0", "./.", "1"):
        unnamed = _prepare(tmp_path, [("A.r1", gt, ".")])
        named = _prepare(tmp_path, [("A.r1", gt, "SVDSS")])
        assert unnamed == named, gt


def test_named_and_unnamed_evidence_in_one_event_is_refused(tmp_path):
    # The unnamed record may come from caller A; counting it as a third
    # caller would give GT 0/1, UC 2, UV 3 instead of ./., 1, 2.
    with pytest.raises(ValueError, match="cannot tell whether they come from the same caller") as error:
        _prepare(tmp_path, [
            ("A.r1", "0/1", "callerA"),
            ("A.r2", "0/0", "callerB"),
            ("A.r3", "0/1", "."),
        ])
    message = str(error.value)
    assert "sampleA.svcf" in message
    assert "A.r1" in message and "A.r3" in message


def _svcf(path, record_id, caller):
    path.write_text(
        "##fileformat=SVCFv1.1\n"
        "##contig=<ID=chr1,length=1000000>\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n"
        f"chr1\t1000\t{record_id}\tN\t<DEL>\t30\tPASS\t"
        "SVTYPE=DEL;END=1500;SVLEN=-500;SUPPORT=10\t"
        "GT:AD:LN:ST:QV:TY:ID:SC:REF:ALT:CO\t"
        f"0/1:5,10:500:+-:30:DEL:{record_id}:{caller}:N:<DEL>:chr1_1000-chr1_1500\n"
    )


def _merge(tmp_path, out_name, callers):
    inputs = []
    for index, caller in enumerate(callers):
        path = tmp_path / out_name / f"s{index}.svcf"
        path.parent.mkdir(exist_ok=True)
        _svcf(path, f"s{index}.1", caller)
        inputs += ["-i", str(path)]
    output = tmp_path / out_name / "merged.svcf"
    result = CliRunner().invoke(
        app,
        ["merge", *inputs, "-o", str(output), "--mode", "sample", "--union"],
    )
    assert result.exit_code == 0, result.output
    return [line for line in output.read_text().splitlines() if not line.startswith("##")]


def test_sample_merge_cli_accepts_single_caller_input_without_caller_name(tmp_path):
    # Different callers in different samples, one of them unnamed.
    unnamed = _merge(tmp_path, "unnamed", ["cuteSV-2.1.3", "."])
    named = _merge(tmp_path, "named", ["cuteSV-2.1.3", "SVDSS"])

    assert len(unnamed) == 2
    assert unnamed == named
