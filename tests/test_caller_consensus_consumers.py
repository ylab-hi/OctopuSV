from pathlib import Path
from types import SimpleNamespace

import pytest

from octopusv.formatter.svcf_to_vcf_converter import SVCFtoVCFConverter
from octopusv.stater.stat_analyzers import GenotypeAnalyzer
from octopusv.stater.stat_reader import SVRecord, read_records
from octopusv.stater.sv_stater import SVStater
from octopusv.utils.caller_consensus import resolve_caller_svcf_consensus


FORMAT = "GT:AD:LN:ST:QV:TY:ID:SC:REF:ALT:CO"


def _segment(gt, ad, ln, record_id, caller):
    return (
        f"{gt}:{ad}:{ln}:.:60:DEL:{record_id}:{caller}:"
        f"N:<DEL>:chr1_100-chr1_200"
    )


def _record(sample_cols, sources="A,B", svlen="100"):
    fields = [
        "chr1",
        "100",
        "event1",
        "N",
        "<DEL>",
        "60",
        "PASS",
        f"SVTYPE=DEL;END=200;SVLEN={svlen};SOURCES={sources}",
        FORMAT,
        *sample_cols,
    ]
    return SVRecord(fields)


def test_caller_consensus_rejects_ambiguous_source_evidence_binding():
    segments = [
        _segment("0/1", "10,8", 100, "a", "A"),
        _segment("1/1", "0,20", 100, "b", "B"),
    ]

    with pytest.raises(ValueError, match="SOURCES/evidence count mismatch"):
        resolve_caller_svcf_consensus(FORMAT, segments, {"SOURCES": "A"})


def test_converter_preserves_single_evidence_without_synthesis():
    converter = object.__new__(SVCFtoVCFConverter)
    segment = _segment("0/1", "10,8", 100, "a", "A")
    event = SimpleNamespace(
        format=FORMAT,
        info={"SOURCES": "A", "SVLEN": "100"},
        sv_type="DEL",
        pos=100,
        end_pos=200,
    )

    assert converter._collapse_caller_blocks(event, [segment]) == "0/1:10,8:18:100"


def test_converter_multi_evidence_uses_consensus_and_drops_noncomposable_ad():
    converter = object.__new__(SVCFtoVCFConverter)
    segments = [
        _segment("0/1", "10,8", 100, "a", "A"),
        _segment("1/1", "0,20", 100, "b", "B"),
    ]
    event = SimpleNamespace(
        format=FORMAT,
        info={"SOURCES": "A,B", "SVLEN": "100"},
        sv_type="DEL",
        pos=100,
        end_pos=200,
    )

    assert converter._collapse_caller_blocks(event, segments) == "1/.:.,.:.:2:2:100"


def test_same_caller_multiple_evidence_is_one_state_not_multiple_votes():
    converter = object.__new__(SVCFtoVCFConverter)
    segments = [
        _segment("0/1", "10,8", 100, "a1", "A"),
        _segment("1/1", "0,20", 100, "a2", "A"),
        _segment("0/1", "12,9", 100, "b", "B"),
    ]
    event = SimpleNamespace(
        format=FORMAT,
        info={"SOURCES": "A,A,B", "SVLEN": "100"},
        sv_type="DEL",
        pos=100,
        end_pos=200,
    )

    assert converter._collapse_caller_blocks(event, segments) == "1/.:.,.:.:2:2:100"



def test_converter_caller_consensus_is_independent_of_evidence_order():
    converter = object.__new__(SVCFtoVCFConverter)
    a = _segment("0/1", "10,8", 100, "a", "A")
    b = _segment("1/1", "0,20", 100, "b", "B")

    event_ab = SimpleNamespace(
        format=FORMAT,
        info={"SOURCES": "A,B", "SVLEN": "100"},
        sv_type="DEL",
        pos=100,
        end_pos=200,
    )
    event_ba = SimpleNamespace(
        format=FORMAT,
        info={"SOURCES": "B,A", "SVLEN": "100"},
        sv_type="DEL",
        pos=100,
        end_pos=200,
    )

    assert converter._collapse_caller_blocks(event_ab, [a, b]) == "1/.:.,.:.:2:2:100"
    assert converter._collapse_caller_blocks(event_ba, [b, a]) == "1/.:.,.:.:2:2:100"

def test_stat_and_converter_use_same_multi_caller_consensus():
    segments = [
        _segment("0/1", "10,8", 100, "a", "A"),
        _segment("1/1", "0,20", 100, "b", "B"),
    ]
    record = _record(segments, sources="A,B", svlen="100")
    event = SimpleNamespace(
        format=FORMAT,
        info={"SOURCES": "A,B", "SVLEN": "100"},
        sv_type="DEL",
        pos=100,
        end_pos=200,
    )

    converter = object.__new__(SVCFtoVCFConverter)
    converted_gt = converter._collapse_caller_blocks(event, segments).split(":", 1)[0]

    stats = GenotypeAnalyzer([record], ["Sample"], mode="caller").analyze()

    assert converted_gt == "1/."
    assert stats["overall"] == {"1/.": 1}


def test_single_sample_multi_marker_is_stat_sample_mode(tmp_path):
    path = Path(tmp_path) / "single_sample_multi.svcf"
    path.write_text(
        "##fileformat=VCFv4.2\n"
        "##OctopuSV_mode=multi\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n"
        "chr1\t100\te1\tN\t<DEL>\t60\tPASS\t"
        "SVTYPE=DEL;END=200;SVLEN=100\t"
        "GT:AD:UC:UV:LN:ST:QV:TY:ID:SC:REF:ALT:CO\t"
        "1/.:.,.:1:2:100:.:60:DEL:e1:OctopuSV:N:<DEL>:chr1_100-chr1_200\n"
    )

    records, sample_names, mode = read_records(path)
    stats = GenotypeAnalyzer(records, sample_names, mode=mode).analyze()

    assert mode == "sample"
    assert stats["mode"] == "sample"
    assert stats["per_sample"] == {"S1": {"1/.": 1}}
    assert stats["overall"] == {"1/.": 1}


def test_svstater_honors_single_sample_multi_marker(tmp_path):
    path = Path(tmp_path) / "single_sample_multi_stater.svcf"
    path.write_text(
        "##fileformat=VCFv4.2\n"
        "##OctopuSV_mode=multi\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n"
        "chr1\t100\te1\tN\t<DEL>\t60\tPASS\t"
        "SVTYPE=DEL;END=200;SVLEN=100\t"
        "GT:AD:UC:UV:LN:ST:QV:TY:ID:SC:REF:ALT:CO\t"
        "1/.:.,.:1:2:100:.:60:DEL:e1:OctopuSV:N:<DEL>:chr1_100-chr1_200\n"
    )

    stater = SVStater(str(path), min_size=0, genome="hg38")
    stater.analyze()

    assert stater.svcf_mode == "sample"
    assert stater.stats["genotype"]["mode"] == "sample"
    assert stater.stats["genotype"]["per_sample"] == {"S1": {"1/.": 1}}
