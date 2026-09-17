from __future__ import annotations

from pathlib import Path

from typer.testing import CliRunner

from octopusv.cli.cli import app
from octopusv.sv import SVEvent
from octopusv.utils.svcf_schema import CALLER_FORMAT, SAMPLE_FORMAT
from octopusv.utils.svcf_utils import write_sv_vcf
from octopusv.utils.svcf_validator import SVCFValidator


RUNNER = CliRunner()


def _sample_block(gt: str, record_id: str, pos: int, end: int) -> str:
    return (
        f"{gt}:.,.:1:1:{end - pos}:+-:60:DEL:{record_id}:OctopuSV:"
        f"N:<DEL>:chr1_{pos}-chr1_{end}"
    )


def _write_sample_mode_svcf(path: Path) -> None:
    records = []
    for pos, s1_gt, s2_gt in [
        (100, "0/1", "1/1"),
        (300, "1/1", "0/1"),
    ]:
        end = pos + 100
        info = (
            f"SVTYPE=DEL;END={end};SVLEN=100;CHR2=chr1;SUPPORT=10;"
            "SVMETHOD=OctopuSV;RTID=.;AF=.;STRAND=+-;RNAMES=.;"
            f"SOURCES=S1,S2;SOURCE_IDS=S1.{pos},S2.{pos}"
        )
        records.append(
            "\t".join(
                [
                    "chr1",
                    str(pos),
                    f"event.{pos}",
                    "N",
                    "<DEL>",
                    "60",
                    "PASS",
                    info,
                    SAMPLE_FORMAT,
                    _sample_block(s1_gt, f"S1.{pos}", pos, end),
                    _sample_block(s2_gt, f"S2.{pos}", pos, end),
                ]
            )
        )

    path.write_text(
        "\n".join(
            [
                "##fileformat=VCFv4.2",
                "##SVCFVersion=1.1",
                "##OctopuSV_mode=multi",
                "##contig=<ID=chr1,length=1000000>",
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2",
                *records,
            ]
        )
        + "\n"
    )


def _write_raw_multi_header(path: Path) -> None:
    path.write_text(
        "\n".join(
            [
                "##fileformat=VCFv4.2",
                "##source=testcaller",
                "##contig=<ID=chr1,length=1000000>",
                '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
                '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths">',
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2",
            ]
        )
        + "\n"
    )


def _legacy_multi_correct_output(path: Path, raw_vcf: Path) -> None:
    event = SVEvent(
        "chr1",
        100,
        "caller.DEL.1",
        "N",
        "<DEL>",
        "60",
        "PASS",
        "SVTYPE=DEL;END=200;SVLEN=-100;CHR2=chr1;STRAND=+-;RE=10",
        format="GT:AD",
        sample="0/1:8,4",
        samples=["0/1:8,4", "0/0:12,0"],
    )
    event.source = "testcaller"
    write_sv_vcf(
        ["##contig=<ID=chr1,length=1000000>"],
        [event],
        path,
        str(raw_vcf),
    )


def test_stat_report_renders_sample_mode_genotype_percentages(tmp_path, monkeypatch):
    monkeypatch.setenv("MPLBACKEND", "Agg")
    input_svcf = tmp_path / "sample.svcf"
    output_txt = tmp_path / "stat.txt"
    _write_sample_mode_svcf(input_svcf)

    result = RUNNER.invoke(
        app,
        ["stat", str(input_svcf), "-o", str(output_txt), "--report"],
    )

    assert result.exit_code == 0, result.output
    output_html = output_txt.with_suffix(".html")
    assert output_html.exists()
    html = output_html.read_text()
    assert "S1 Genotypes" in html
    assert "S2 Genotypes" in html
    assert "50.00%" in html


def test_validator_accepts_correct_legacy_multi_with_warning(tmp_path):
    raw_vcf = tmp_path / "multi.vcf"
    output_svcf = tmp_path / "multi.svcf"
    _write_raw_multi_header(raw_vcf)
    _legacy_multi_correct_output(output_svcf, raw_vcf)

    validator = SVCFValidator(str(output_svcf))
    validator.validate()

    assert validator.status() == "PASS_WITH_WARNINGS"
    assert validator.exit_code() == 0
    assert validator.blocking is False
    assert not validator.errors
    assert [issue.code for issue in validator.warnings] == ["W_LEGACY_MULTI_001"]

    cli_result = RUNNER.invoke(app, ["validate-svcf", "-i", str(output_svcf)])
    assert cli_result.exit_code == 0, cli_result.output
    assert "Status: PASS_WITH_WARNINGS" in cli_result.output
    assert "Blocking for Downstream: False" in cli_result.output


def test_validator_keeps_versioned_multi_strict(tmp_path):
    path = tmp_path / "invalid_v11_multi.svcf"
    block = (
        "0/1:5,7:100:+-:60:DEL:source.1:testcaller:N:<DEL>:"
        "chr1_100-chr1_200"
    )
    info = (
        "SVTYPE=DEL;END=200;SVLEN=100;CHR2=chr1;SUPPORT=5;"
        "SVMETHOD=OctopuSV;RTID=.;AF=.;STRAND=+-;RNAMES=."
    )
    path.write_text(
        "\n".join(
            [
                "##fileformat=VCFv4.2",
                "##SVCFVersion=1.1",
                "##OctopuSV_mode=multi",
                "##contig=<ID=chr1,length=1000000>",
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2",
                f"chr1\t100\te1\tN\t<DEL>\t60\tPASS\t{info}\t{CALLER_FORMAT}\t{block}\t{block}",
            ]
        )
        + "\n"
    )

    validator = SVCFValidator(str(path))
    validator.validate()

    assert validator.status() == "FAILED"
    assert validator.blocking is True
    assert "E_FMT_001" in {issue.code for issue in validator.errors}


def test_validator_rejects_mixed_legacy_multi_schemas(tmp_path):
    path = tmp_path / "mixed_legacy_multi.svcf"
    caller_block = (
        "0/1:5,7:100:+-:60:DEL:caller.1:testcaller:N:<DEL>:"
        "chr1_100-chr1_200"
    )
    sample_block = _sample_block("0/1", "sample.1", 300, 400)
    info = (
        "SVTYPE=DEL;END=200;SVLEN=100;CHR2=chr1;SUPPORT=5;"
        "SVMETHOD=OctopuSV;RTID=.;AF=.;STRAND=+-;RNAMES=."
    )
    info2 = info.replace("END=200", "END=400")
    path.write_text(
        "\n".join(
            [
                "##fileformat=VCFv4.2",
                "##OctopuSV_mode=multi",
                "##contig=<ID=chr1,length=1000000>",
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2",
                f"chr1\t100\te1\tN\t<DEL>\t60\tPASS\t{info}\t{CALLER_FORMAT}\t{caller_block}\t{caller_block}",
                f"chr1\t300\te2\tN\t<DEL>\t60\tPASS\t{info2}\t{SAMPLE_FORMAT}\t{sample_block}\t{sample_block}",
            ]
        )
        + "\n"
    )

    validator = SVCFValidator(str(path))
    validator.validate()

    assert validator.status() == "FAILED"
    assert "E_FMT_001" in {issue.code for issue in validator.errors}
