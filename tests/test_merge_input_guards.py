from pathlib import Path
from types import SimpleNamespace

import pytest

from octopusv.cli.merge import (
    _preserve_sample_input_evidence,
    _preflight_merge_inputs,
)


def _write_svcf(path: Path, *, header_samples, evidence_blocks, multi_marker=False):
    path.parent.mkdir(parents=True, exist_ok=True)
    header = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT"
    if header_samples:
        header += "\t" + "\t".join(header_samples)

    fields = [
        "chr1",
        "100",
        "event1",
        "N",
        "<INS>",
        "60",
        "PASS",
        "SVTYPE=INS;END=150;SVLEN=50;CHR2=chr1;SUPPORT=5;SVMETHOD=OctopuSV;RTID=.;AF=.;STRAND=.;RNAMES=.",
        "GT:AD:LN:ST:QV:TY:ID:SC:REF:ALT:CO",
    ]
    fields.extend(evidence_blocks)

    meta = "##fileformat=VCFv4.2\n"
    if multi_marker:
        meta += "##OctopuSV_mode=multi\n"

    path.write_text(
        meta
        + header
        + "\n"
        + "\t".join(fields)
        + "\n"
    )


BLOCK_A = "0/1:5,5:50:.:60:INS:a.1:callerA:N:<INS>:chr1_100-chr1_150"
BLOCK_B = "0/1:4,6:50:.:55:INS:b.1:callerB:N:<INS>:chr1_102-chr1_152"


def test_preflight_rejects_duplicate_realpath(tmp_path):
    real_file = tmp_path / "real.svcf"
    _write_svcf(real_file, header_samples=["SAMPLE"], evidence_blocks=[BLOCK_A])

    link_file = tmp_path / "alias.svcf"
    link_file.symlink_to(real_file)

    with pytest.raises(ValueError, match="same physical file"):
        _preflight_merge_inputs(
            input_files=[real_file, link_file],
            labels=["real", "alias"],
            mode="caller",
        )


def test_preflight_rejects_duplicate_labels_with_actionable_option(tmp_path):
    file_a = tmp_path / "sampleA" / "merged.svcf"
    file_b = tmp_path / "sampleB" / "merged.svcf"
    _write_svcf(file_a, header_samples=["SAMPLE"], evidence_blocks=[BLOCK_A])
    _write_svcf(file_b, header_samples=["SAMPLE"], evidence_blocks=[BLOCK_B])

    with pytest.raises(ValueError) as exc_info:
        _preflight_merge_inputs(
            input_files=[file_a, file_b],
            labels=["merged", "merged"],
            mode="sample",
        )

    message = str(exc_info.value)
    assert str(file_a) in message
    assert str(file_b) in message
    assert "--sample-names" in message


def test_caller_mode_rejects_multi_evidence_input(tmp_path):
    merged = tmp_path / "merged.svcf"
    _write_svcf(
        merged,
        header_samples=["SAMPLE"],
        evidence_blocks=[BLOCK_A, BLOCK_B],
    )

    with pytest.raises(ValueError, match="Caller-mode re-merge"):
        _preflight_merge_inputs(
            input_files=[merged],
            labels=["merged"],
            mode="caller",
        )


def test_sample_mode_allows_caller_merged_single_sample_input(tmp_path):
    """Protect the caller-merge-per-sample -> population-merge workflow."""
    merged = tmp_path / "sample1.svcf"
    _write_svcf(
        merged,
        header_samples=["SAMPLE"],
        evidence_blocks=[BLOCK_A, BLOCK_B],
    )

    _preflight_merge_inputs(
        input_files=[merged],
        labels=["sample1"],
        mode="sample",
    )


def test_sample_mode_rejects_true_multi_sample_input(tmp_path):
    population = tmp_path / "population.svcf"
    _write_svcf(
        population,
        header_samples=["sample1", "sample2"],
        evidence_blocks=[BLOCK_A, BLOCK_B],
    )

    with pytest.raises(ValueError, match="expects one input file per biological sample"):
        _preflight_merge_inputs(
            input_files=[population],
            labels=["population"],
            mode="sample",
        )




def test_merge_rejects_synthesized_sample_multi_input_even_with_one_sample(tmp_path):
    """A one-column synthesized sample SVCF is still not caller evidence."""
    population_subset = tmp_path / "population_subset.svcf"
    _write_svcf(
        population_subset,
        header_samples=["sample1"],
        evidence_blocks=[BLOCK_A],
        multi_marker=True,
    )

    with pytest.raises(ValueError, match="synthesized sample/multi SVCF"):
        _preflight_merge_inputs(
            input_files=[population_subset],
            labels=["sample1"],
            mode="sample",
        )


def test_caller_merge_also_rejects_synthesized_sample_multi_input(tmp_path):
    """A synthesized sample call must not masquerade as caller evidence."""
    population_subset = tmp_path / "population_subset.svcf"
    _write_svcf(
        population_subset,
        header_samples=["sample1"],
        evidence_blocks=[BLOCK_A],
        multi_marker=True,
    )

    with pytest.raises(ValueError, match="synthesized sample/multi SVCF"):
        _preflight_merge_inputs(
            input_files=[population_subset],
            labels=["sample1"],
            mode="caller",
        )


def test_preserve_sample_input_evidence_keeps_full_payload():
    sample = {"ID": "a.1", "GT": "0/1"}
    event = SimpleNamespace(
        sample=sample,
        format="GT:AD:LN:ST:QV:TY:ID:SC:REF:ALT:CO",
        info={
            "SOURCES": "cuteSV,svim,svim",
            "SOURCE_IDS": "a.1,b.1,b.2",
        },
        raw_sample_columns=[BLOCK_A, BLOCK_B, BLOCK_B],
    )

    _preserve_sample_input_evidence([event])

    assert sample["ID"] == "a.1"
    assert "_octopusv_collapsed_evidence_count" not in sample

    payload = sample["_octopusv_evidence_payload"]
    assert payload["format"] == "GT:AD:LN:ST:QV:TY:ID:SC:REF:ALT:CO"
    assert payload["blocks"] == (BLOCK_A, BLOCK_B, BLOCK_B)
    assert payload["sources"] == "cuteSV,svim,svim"
    assert payload["source_ids"] == "a.1,b.1,b.2"
    assert payload["record"]["id"] == "."


def test_preserve_sample_input_evidence_keeps_single_evidence_too():
    sample = {"ID": "a.1", "GT": "0/1"}
    event = SimpleNamespace(
        sample=sample,
        format="GT:AD:LN:ST:QV:TY:ID:SC:REF:ALT:CO",
        info={},
        raw_sample_columns=[BLOCK_A],
    )

    _preserve_sample_input_evidence([event])

    assert "_octopusv_collapsed_evidence_count" not in sample
    assert sample["_octopusv_evidence_payload"]["blocks"] == (BLOCK_A,)


def test_caller_mode_preflight_scans_beyond_first_data_record(tmp_path):
    """A later multi-evidence record must not escape the caller-mode guard."""
    path = tmp_path / "late_multi.svcf"
    header = (
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n"
    )
    one = (
        "chr1\t100\tevent{idx}\tN\t<INS>\t60\tPASS\t"
        "SVTYPE=INS;END=150;SVLEN=50;CHR2=chr1;SUPPORT=5;"
        "SVMETHOD=OctopuSV;RTID=.;AF=.;STRAND=.;RNAMES=.\t"
        "GT:AD:LN:ST:QV:TY:ID:SC:REF:ALT:CO\t" + BLOCK_A + "\n"
    )
    lines = [header]
    for idx in range(1, 50):
        lines.append(one.format(idx=idx))
    lines.append(
        "chr1\t200\tlate\tN\t<INS>\t60\tPASS\t"
        "SVTYPE=INS;END=250;SVLEN=50;CHR2=chr1;SUPPORT=5;"
        "SVMETHOD=OctopuSV;RTID=.;AF=.;STRAND=.;RNAMES=.\t"
        "GT:AD:LN:ST:QV:TY:ID:SC:REF:ALT:CO\t"
        + BLOCK_A + "\t" + BLOCK_B + "\n"
    )
    path.write_text("".join(lines))

    with pytest.raises(ValueError, match="Caller-mode re-merge"):
        _preflight_merge_inputs(
            input_files=[path],
            labels=["late_multi"],
            mode="caller",
        )


def test_preflight_rejects_missing_chrom_header(tmp_path):
    path = tmp_path / "missing_header.svcf"
    path.write_text(
        "##fileformat=VCFv4.2\n"
        "chr1\t100\tevent1\tN\t<INS>\t60\tPASS\t"
        "SVTYPE=INS;END=150;SVLEN=50;CHR2=chr1;SUPPORT=5;"
        "SVMETHOD=OctopuSV;RTID=.;AF=.;STRAND=.;RNAMES=.\t"
        "GT:AD:LN:ST:QV:TY:ID:SC:REF:ALT:CO\t"
        + BLOCK_A + "\n"
    )

    with pytest.raises(ValueError, match="missing a #CHROM header"):
        _preflight_merge_inputs(
            input_files=[path],
            labels=["missing_header"],
            mode="caller",
        )


def test_preflight_rejects_unknown_specific_before_shape_scan(tmp_path, monkeypatch):
    input_file = tmp_path / "input.svcf"
    unknown = tmp_path / "unknown.svcf"
    _write_svcf(input_file, header_samples=["SAMPLE"], evidence_blocks=[BLOCK_A])

    import octopusv.cli.merge as merge_module

    calls = []
    original = merge_module._preflight_svcf_shape

    def counted_shape(path, mode):
        calls.append(str(path))
        return original(path, mode)

    monkeypatch.setattr(merge_module, "_preflight_svcf_shape", counted_shape)

    with pytest.raises(ValueError, match="not one of the merge inputs"):
        _preflight_merge_inputs(
            input_files=[input_file],
            labels=["input"],
            mode="caller",
            specific=[unknown],
        )

    assert calls == []


def test_preflight_rejects_expression_basename_collision_before_shape_scan(tmp_path, monkeypatch):
    source_a = tmp_path / "d1" / "shared.svcf"
    source_b = tmp_path / "d2" / "shared.svcf"
    _write_svcf(source_a, header_samples=["SAMPLE"], evidence_blocks=[BLOCK_A])
    _write_svcf(source_b, header_samples=["SAMPLE"], evidence_blocks=[BLOCK_B])

    import octopusv.cli.merge as merge_module

    calls = []
    original = merge_module._preflight_svcf_shape

    def counted_shape(path, mode):
        calls.append(str(path))
        return original(path, mode)

    monkeypatch.setattr(merge_module, "_preflight_svcf_shape", counted_shape)

    with pytest.raises(ValueError, match="share the basename"):
        _preflight_merge_inputs(
            input_files=[source_a, source_b],
            labels=["a", "b"],
            mode="caller",
            expression="shared.svcf",
        )

    assert calls == []


def test_preflight_rejects_expression_identifier_collision_before_shape_scan(tmp_path, monkeypatch):
    source_a = tmp_path / "a-b.svcf"
    source_b = tmp_path / "a_b.svcf"
    _write_svcf(source_a, header_samples=["SAMPLE"], evidence_blocks=[BLOCK_A])
    _write_svcf(source_b, header_samples=["SAMPLE"], evidence_blocks=[BLOCK_B])

    import octopusv.cli.merge as merge_module

    calls = []
    original = merge_module._preflight_svcf_shape

    def counted_shape(path, mode):
        calls.append(str(path))
        return original(path, mode)

    monkeypatch.setattr(merge_module, "_preflight_svcf_shape", counted_shape)

    with pytest.raises(ValueError, match="same expression identifier"):
        _preflight_merge_inputs(
            input_files=[source_a, source_b],
            labels=["a", "b"],
            mode="caller",
            expression="a-b.svcf OR a_b.svcf",
        )

    assert calls == []
