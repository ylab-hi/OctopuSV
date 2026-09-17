from __future__ import annotations

import json
from pathlib import Path

import pytest
from typer.testing import CliRunner

pytest.importorskip("natsort")

from octopusv.bencher.sv_bencher import SVBencher
from octopusv.cli.cli import app
from octopusv.subset.svcf_subset import SVCFSubset, SubsetConfig
from octopusv.utils.svcf_schema import CALLER_FORMAT, SAMPLE_FORMAT
from octopusv.utils.svcf_validator import SVCFValidator


runner = CliRunner()


def _caller_block(record_id: str, source: str, *, pos: int = 100, end: int = 200, gt: str = "0/1") -> str:
    return (
        f"{gt}:5,5:{end-pos}:.:60:DEL:{record_id}:{source}:N:<DEL>:"
        f"chr1_{pos}-chr1_{end}"
    )


def _sample_block(record_id: str, source: str, *, pos: int = 100, end: int = 200, gt: str = "0/1") -> str:
    return (
        f"{gt}:5,5:1:1:{end-pos}:.:60:DEL:{record_id}:{source}:N:<DEL>:"
        f"chr1_{pos}-chr1_{end}"
    )


def _info(*, end: str = "200", sources: str | None = None, source_ids: str | None = None) -> str:
    items = [
        "SVTYPE=DEL",
        f"END={end}",
        "SVLEN=100",
        "CHR2=chr1",
        "SUPPORT=5",
        "SVMETHOD=OctopuSV",
        "RTID=.",
        "AF=.",
        "STRAND=.",
        "RNAMES=.",
    ]
    if sources is not None:
        items.append(f"SOURCES={sources}")
    if source_ids is not None:
        items.append(f"SOURCE_IDS={source_ids}")
    return ";".join(items)


def _write_caller_svcf(
    path: Path,
    *,
    source: str = "PBSV",
    end: str = "200",
    versioned: bool = True,
    include_custom_meta: bool = False,
    include_sources: bool = True,
) -> None:
    meta = ["##fileformat=VCFv4.2"]
    if versioned:
        meta += ["##SVCFVersion=1.1", "##OctopuSV_mode=caller"]
    meta += [
        "##contig=<ID=chr1,length=1000000>",
        '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="SV type">',
        '##INFO=<ID=END,Number=1,Type=Integer,Description="End">',
        '##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Length">',
        '##INFO=<ID=CHR2,Number=1,Type=String,Description="Mate contig">',
        '##INFO=<ID=SUPPORT,Number=1,Type=Integer,Description="Support">',
        '##INFO=<ID=SVMETHOD,Number=1,Type=String,Description="Method">',
        '##INFO=<ID=RTID,Number=1,Type=String,Description="Related ID">',
        '##INFO=<ID=AF,Number=1,Type=String,Description="AF">',
        '##INFO=<ID=STRAND,Number=1,Type=String,Description="Strand">',
        '##INFO=<ID=RNAMES,Number=1,Type=String,Description="Reads">',
        '##INFO=<ID=SOURCES,Number=.,Type=String,Description="Sources">',
        '##INFO=<ID=SOURCE_IDS,Number=.,Type=String,Description="Source IDs">',
    ]
    if include_custom_meta:
        meta += [
            '##FILTER=<ID=LowQual,Description="Low quality">',
            '##ALT=<ID=DEL,Description="Deletion">',
            '##INFO=<ID=CUSTOM,Number=1,Type=String,Description="Custom annotation">',
            "##reference=GRCh38",
            "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
        ]
    record_id = "r1"
    info = _info(
        end=end,
        sources=source if include_sources else None,
        source_ids=record_id if include_sources else None,
    )
    if include_custom_meta:
        info += ";CUSTOM=abc"
    block = _caller_block(record_id, source, end=int(end) if end.isdigit() else 200)
    path.write_text(
        "\n".join(
            meta
            + [
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE",
                f"chr1\t100\t{record_id}\tN\t<DEL>\t60\tPASS\t{info}\t{CALLER_FORMAT}\t{block}",
            ]
        )
        + "\n"
    )


def _write_sample_svcf(path: Path, *, versioned: bool, truncate_record: bool) -> None:
    meta = ["##fileformat=VCFv4.2"]
    if versioned:
        meta += ["##SVCFVersion=1.1", "##OctopuSV_mode=multi"]
    row_blocks = [_sample_block("r1", "S1"), _sample_block("r1", "S2")]
    if truncate_record:
        row_blocks = row_blocks[:1]
    path.write_text(
        "\n".join(
            meta
            + [
                "##contig=<ID=chr1,length=1000000>",
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2",
                "chr1\t100\tr1\tN\t<DEL>\t60\tPASS\t"
                + _info(sources="S1,S2", source_ids="r1,r1")
                + f"\t{SAMPLE_FORMAT}\t"
                + "\t".join(row_blocks),
            ]
        )
        + "\n"
    )


def _data_rows(path: Path) -> list[str]:
    return [line for line in path.read_text().splitlines() if line and not line.startswith("#")]


def test_subset_case_only_caller_typo_fails_with_suggestion(tmp_path):
    source = tmp_path / "input.svcf"
    output = tmp_path / "out.svcf"
    _write_caller_svcf(source, source="PBSV")

    result = runner.invoke(
        app,
        ["subset", "-i", str(source), "-o", str(output), "--caller", "pbsv"],
    )

    assert result.exit_code != 0
    assert "Did you mean 'PBSV'" in result.output
    assert not output.exists()


def test_subset_unobserved_caller_warns_and_returns_empty_success(tmp_path):
    source = tmp_path / "input.svcf"
    output = tmp_path / "out.svcf"
    _write_caller_svcf(source, source="A")

    result = runner.invoke(
        app,
        ["subset", "-i", str(source), "-o", str(output), "--caller", "ZZZ"],
    )

    assert result.exit_code == 0, result.output
    assert "Requested caller 'ZZZ' was not observed" in result.output
    assert "Observed callers: A" in result.output
    assert _data_rows(output) == []
    validator = SVCFValidator(str(output))
    validator.validate()
    assert not validator.blocking, validator.to_summary(max_issues=None)


def test_subset_single_evidence_sc_is_an_explicit_caller_identity(tmp_path):
    source = tmp_path / "single.svcf"
    output = tmp_path / "out.svcf"
    _write_caller_svcf(source, source="PBSV", include_sources=False)

    result = runner.invoke(
        app,
        ["subset", "-i", str(source), "-o", str(output), "--caller", "PBSV"],
    )

    assert result.exit_code == 0, result.output
    rows = _data_rows(output)
    assert len(rows) == 1
    assert "SOURCES=PBSV" in rows[0]


def test_subset_multiple_callers_keeps_observed_and_warns_for_unobserved(tmp_path):
    source = tmp_path / "input.svcf"
    output = tmp_path / "out.svcf"
    _write_caller_svcf(source, source="A")

    result = runner.invoke(
        app,
        ["subset", "-i", str(source), "-o", str(output), "--caller", "A,ZZZ"],
    )

    assert result.exit_code == 0, result.output
    assert "Requested caller 'ZZZ' was not observed" in result.output
    assert "Observed callers: A" in result.output
    assert len(_data_rows(output)) == 1
    assert "SOURCES=A" in _data_rows(output)[0]


def test_subset_versioned_v11_rejects_sample_column_mismatch_before_processing(tmp_path):
    source = tmp_path / "bad.svcf"
    output = tmp_path / "out.svcf"
    _write_sample_svcf(source, versioned=True, truncate_record=True)

    result = runner.invoke(
        app,
        ["subset", "-i", str(source), "-o", str(output), "--sample", "S1"],
    )

    assert result.exit_code != 0
    assert "E_COL_001" in result.output
    assert not output.exists()


def test_subset_versioned_v11_rejects_caller_source_evidence_mismatch(tmp_path):
    source = tmp_path / "bad_caller.svcf"
    output = tmp_path / "out.svcf"
    block = _caller_block("r1", "A")
    source.write_text(
        "##fileformat=VCFv4.2\n"
        "##SVCFVersion=1.1\n"
        "##OctopuSV_mode=caller\n"
        "##contig=<ID=chr1,length=1000000>\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n"
        "chr1\t100\tr1\tN\t<DEL>\t60\tPASS\t"
        + _info(sources="A,B", source_ids="r1,r2")
        + f"\t{CALLER_FORMAT}\t{block}\n"
    )

    result = runner.invoke(
        app,
        ["subset", "-i", str(source), "-o", str(output), "--caller", "A"],
    )

    assert result.exit_code != 0
    assert "E_COL_002" in result.output
    assert not output.exists()


def test_subset_legacy_column_mismatch_keeps_compatibility_behavior(tmp_path):
    source = tmp_path / "legacy.svcf"
    output = tmp_path / "out.svcf"
    _write_sample_svcf(source, versioned=False, truncate_record=True)

    result = runner.invoke(
        app,
        ["subset", "-i", str(source), "-o", str(output), "--sample", "S1"],
    )

    assert result.exit_code == 0, result.output
    assert len(_data_rows(output)) == 1


def test_benchmark_rejects_bad_versioned_v11_before_metrics(tmp_path):
    truth = tmp_path / "truth.svcf"
    call = tmp_path / "call.svcf"
    out = tmp_path / "bench"
    _write_caller_svcf(truth, end="oops", versioned=True)
    _write_caller_svcf(call, end="200", versioned=True)

    with pytest.raises(ValueError, match=r"E_END_002"):
        SVBencher(truth, call, out).run_benchmark()

    assert not (out / "summary.json").exists()


def test_benchmark_legacy_bad_end_keeps_legacy_compatibility(tmp_path):
    truth = tmp_path / "truth.svcf"
    call = tmp_path / "call.svcf"
    out = tmp_path / "bench"
    _write_caller_svcf(truth, end="oops", versioned=False)
    _write_caller_svcf(call, end="200", versioned=False)

    SVBencher(truth, call, out).run_benchmark()
    assert (out / "summary.json").exists()


def test_benchmark_vcf_outputs_are_self_describing_and_inherit_safe_header(tmp_path):
    truth = tmp_path / "truth.svcf"
    call = tmp_path / "call.svcf"
    out = tmp_path / "bench"
    _write_caller_svcf(truth, include_custom_meta=True)
    _write_caller_svcf(call, include_custom_meta=True)

    SVBencher(truth, call, out).run_benchmark()

    for name in ["tp-base.vcf", "tp-call.vcf", "fp.vcf", "fn.vcf"]:
        text = (out / name).read_text()
        assert text.startswith("##fileformat=VCFv4.2\n")
        assert "##contig=<ID=chr1,length=1000000>" in text
        assert "##INFO=<ID=CUSTOM," in text
        assert "##FILTER=<ID=LowQual," in text
        assert "##ALT=<ID=DEL," in text
        assert "##reference=GRCh38" in text
        assert "##SVCFVersion=" not in text
        assert "##OctopuSV_mode=" not in text
        assert "##FORMAT=" not in text
        declared = {
            line.split("##INFO=<ID=", 1)[1].split(",", 1)[0]
            for line in text.splitlines()
            if line.startswith("##INFO=<ID=")
        }
        for row in _data_rows(out / name):
            used = {
                item.split("=", 1)[0]
                for item in row.split("\t")[7].split(";")
                if item
            }
            assert used <= declared

    metrics = json.loads((out / "summary.json").read_text())
    assert metrics["TP"] == 1
    assert metrics["FP"] == 0
    assert metrics["FN"] == 0
    assert metrics["F1"] == 1.0


def _write_empty_svcf(path: Path) -> None:
    path.write_text(
        "##fileformat=VCFv4.2\n"
        "##SVCFVersion=1.1\n"
        "##OctopuSV_mode=caller\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n"
    )


@pytest.mark.parametrize("command,suffix", [("svcf2bed", "bed"), ("svcf2bedpe", "bedpe")])
def test_empty_svcf_conversion_is_successful(command, suffix, tmp_path):
    source = tmp_path / "empty.svcf"
    output = tmp_path / f"out.{suffix}"
    _write_empty_svcf(source)

    result = runner.invoke(app, [command, "-i", str(source), "-o", str(output)])

    assert result.exit_code == 0, result.output
    assert output.exists()
    assert "Error occurred:" not in result.output
    if command == "svcf2bedpe":
        assert output.read_text().startswith("#chrom1\tstart1\tend1\tchrom2")
    else:
        assert output.read_text().startswith("track name=SVs")


def test_source_label_colon_is_rejected_by_validator_and_producers(tmp_path):
    # Explicit versioned SOURCES atom.
    source = tmp_path / "bad_source.svcf"
    _write_caller_svcf(source, source="A:B")
    validator = SVCFValidator(str(source))
    validator.validate()
    assert "E_SRC_005" in {issue.code for issue in validator.issues}

    # correct must fail before serializing an ambiguous SC block.
    raw = tmp_path / "raw.vcf"
    corrected = tmp_path / "corrected.svcf"
    raw.write_text(
        "##fileformat=VCFv4.2\n"
        "##source=A:B\n"
        "##contig=<ID=chr1,length=1000000>\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n"
        "chr1\t100\tr1\tN\t<DEL>\t60\tPASS\tSVTYPE=DEL;END=200;SVLEN=-100\tGT\t0/1\n"
    )
    result = runner.invoke(app, ["correct", "-i", str(raw), "-o", str(corrected)])
    assert result.exit_code != 0
    assert "A:B" in result.output
    assert not corrected.exists()


def test_merge_caller_name_with_colon_is_rejected(tmp_path):
    source = tmp_path / "input.svcf"
    output = tmp_path / "out.svcf"
    _write_caller_svcf(source, source="A")

    result = runner.invoke(
        app,
        [
            "merge",
            "-i",
            str(source),
            "-o",
            str(output),
            "--mode",
            "caller",
            "--caller-names",
            "A:B",
            "--union",
        ],
    )
    assert result.exit_code != 0
    assert "A:B" in result.output
    assert not output.exists()
