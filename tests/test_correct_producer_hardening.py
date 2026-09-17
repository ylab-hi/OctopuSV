from __future__ import annotations

from pathlib import Path

import pytest
from typer.testing import CliRunner

pytest.importorskip("natsort")

from octopusv.cli.cli import app
from octopusv.utils.svcf_validator import SVCFValidator


runner = CliRunner()


def _header(*, fmt_lines: list[str] | None = None) -> str:
    lines = [
        "##fileformat=VCFv4.2",
        "##source=TestCaller",
        "##contig=<ID=chr1,length=1000000>",
        "##contig=<ID=chr5,length=1000000>",
    ]
    if fmt_lines is None:
        fmt_lines = [
            '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
        ]
    lines.extend(fmt_lines)
    lines.append("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1")
    return "\n".join(lines) + "\n"


def _record(
    *,
    chrom: str = "chr1",
    pos: int = 100,
    record_id: str = "r1",
    alt: str = "<DEL>",
    info: str = "SVTYPE=DEL;END=150;SVLEN=-50",
    fmt: str = "GT",
    sample: str = "0/1",
) -> str:
    return (
        f"{chrom}\t{pos}\t{record_id}\tN\t{alt}\t60\tPASS\t"
        f"{info}\t{fmt}\t{sample}\n"
    )


def _invoke_correct(
    input_vcf: Path,
    output_svcf: Path,
    *extra_args: str,
):
    return runner.invoke(
        app,
        [
            "correct",
            "-i",
            str(input_vcf),
            "-o",
            str(output_svcf),
            *extra_args,
        ],
    )


def _data_lines(path: Path) -> list[str]:
    return [
        line
        for line in path.read_text().splitlines()
        if line and not line.startswith("#")
    ]


def _info_dict(record_line: str) -> dict[str, str]:
    info = record_line.split("\t")[7]
    result: dict[str, str] = {}
    for item in info.split(";"):
        if "=" in item:
            key, value = item.split("=", 1)
            result[key] = value
    return result


def _assert_valid_svcf(path: Path) -> None:
    validator = SVCFValidator(str(path))
    validator.validate()
    assert validator.errors == [], validator.to_summary(max_issues=None)


def test_correct_rejects_illegal_final_svtype_without_output(tmp_path):
    input_vcf = tmp_path / "cpx.vcf"
    output_svcf = tmp_path / "out.svcf"
    input_vcf.write_text(
        _header()
        + _record(
            alt="<CPX>",
            info="SVTYPE=CPX;END=150;SVLEN=50",
        )
    )

    result = _invoke_correct(input_vcf, output_svcf)

    assert result.exit_code != 0
    assert "CPX" in result.output
    assert "SVTYPE" in result.output
    assert not output_svcf.exists()


def test_correct_derives_missing_span_end_from_svlen_and_validates(tmp_path):
    input_vcf = tmp_path / "missing_end.vcf"
    output_svcf = tmp_path / "out.svcf"
    input_vcf.write_text(
        _header()
        + _record(
            pos=100,
            info="SVTYPE=DEL;SVLEN=-50",
        )
    )

    result = _invoke_correct(input_vcf, output_svcf)

    assert result.exit_code == 0, result.output
    record = _data_lines(output_svcf)[0]
    info = _info_dict(record)
    assert info["SVTYPE"] == "DEL"
    assert info["END"] == "150"
    assert info["SVLEN"] == "50"
    _assert_valid_svcf(output_svcf)


def test_correct_rejects_malformed_explicit_span_end_instead_of_guessing(tmp_path):
    input_vcf = tmp_path / "bad_end.vcf"
    output_svcf = tmp_path / "out.svcf"
    input_vcf.write_text(
        _header()
        + _record(
            info="SVTYPE=DEL;END=oops;SVLEN=-50",
        )
    )

    result = _invoke_correct(input_vcf, output_svcf)

    assert result.exit_code != 0
    assert "END" in result.output
    assert "oops" in result.output
    assert not output_svcf.exists()


def test_correct_preserves_symbolic_tra_with_known_breakpoints(tmp_path):
    input_vcf = tmp_path / "symbolic_tra.vcf"
    output_svcf = tmp_path / "out.svcf"
    input_vcf.write_text(
        _header()
        + _record(
            pos=100,
            record_id="tra1",
            alt="<TRA>",
            info="SVTYPE=TRA;CHR2=chr5;END=5000",
        )
    )

    result = _invoke_correct(input_vcf, output_svcf)

    assert result.exit_code == 0, result.output
    fields = _data_lines(output_svcf)[0].split("\t")
    info = _info_dict("\t".join(fields))
    assert fields[4] == "<TRA>"
    assert info["SVTYPE"] == "TRA"
    assert info["CHR2"] == "chr5"
    assert info["END"] == "5000"
    assert info["SVLEN"] == "."
    _assert_valid_svcf(output_svcf)


def test_correct_recovers_bracket_tra_coordinates_from_alt(tmp_path):
    input_vcf = tmp_path / "bracket_tra_missing_coords.vcf"
    output_svcf = tmp_path / "out.svcf"
    input_vcf.write_text(
        _header()
        + _record(
            pos=100,
            record_id="tra_bracket",
            alt="N]chr5:5000]",
            info="SVTYPE=TRA",
        )
    )

    result = _invoke_correct(input_vcf, output_svcf)

    assert result.exit_code == 0, result.output
    record = _data_lines(output_svcf)[0]
    info = _info_dict(record)
    assert info["SVTYPE"] == "TRA"
    assert info["CHR2"] == "chr5"
    assert info["END"] == "5000"
    assert info["SVLEN"] == "."
    _assert_valid_svcf(output_svcf)


def test_correct_does_not_invent_chr2_for_symbolic_tra(tmp_path):
    input_vcf = tmp_path / "symbolic_tra_missing_chr2.vcf"
    output_svcf = tmp_path / "out.svcf"
    input_vcf.write_text(
        _header()
        + _record(
            pos=100,
            record_id="tra_missing_chr2",
            alt="<TRA>",
            info="SVTYPE=TRA;END=5000",
        )
    )

    result = _invoke_correct(input_vcf, output_svcf)

    assert result.exit_code != 0
    assert "CHR2" in result.output
    assert "orientation may remain unknown" in result.output
    assert not output_svcf.exists()


def test_correct_rejects_clearly_structural_record_missing_svtype(tmp_path):
    input_vcf = tmp_path / "missing_svtype.vcf"
    output_svcf = tmp_path / "out.svcf"
    input_vcf.write_text(
        _header()
        + _record(
            alt="<DEL>",
            info="END=150;SVLEN=-50",
        )
    )

    result = _invoke_correct(input_vcf, output_svcf)

    assert result.exit_code != 0
    assert "SVTYPE" in result.output
    assert "<DEL>" in result.output
    assert not output_svcf.exists()


def test_correct_keeps_reference_block_compatibility_when_real_sv_exists(tmp_path):
    input_vcf = tmp_path / "mixed_reference_and_sv.vcf"
    output_svcf = tmp_path / "out.svcf"
    reference_block = _record(
        pos=1,
        record_id="refblock",
        alt="<NON_REF>",
        info="END=99",
        sample="0/0",
    )
    real_sv = _record(
        pos=100,
        record_id="del1",
        alt="<DEL>",
        info="SVTYPE=DEL;END=150;SVLEN=-50",
    )
    input_vcf.write_text(_header() + reference_block + real_sv)

    result = _invoke_correct(input_vcf, output_svcf)

    assert result.exit_code == 0, result.output
    records = _data_lines(output_svcf)
    assert len(records) == 1
    assert records[0].split("\t")[2] == "del1"
    assert "Skipped 1 non-SV/reference record" in result.output
    _assert_valid_svcf(output_svcf)


def test_correct_refuses_silent_empty_output_when_input_has_no_sv_records(tmp_path):
    input_vcf = tmp_path / "reference_only.vcf"
    output_svcf = tmp_path / "out.svcf"
    input_vcf.write_text(
        _header()
        + _record(
            pos=1,
            record_id="refblock",
            alt="<NON_REF>",
            info="END=100",
            sample="0/0",
        )
    )

    result = _invoke_correct(input_vcf, output_svcf)

    assert result.exit_code != 0
    assert "no structural-variant records" in result.output.lower()
    assert not output_svcf.exists()


def test_correct_validation_backstop_preserves_existing_output(tmp_path, monkeypatch):
    input_vcf = tmp_path / "input.vcf"
    output_svcf = tmp_path / "out.svcf"
    input_vcf.write_text(_header() + _record())
    output_svcf.write_text("ORIGINAL\n")

    def invalid_but_successful_writer(
        contig_lines,
        events,
        output_file,
        input_vcf_file=None,
        extra_meta_lines=None,
    ):
        Path(output_file).write_text(
            "##fileformat=VCFv4.2\n"
            "##SVCFVersion=1.1\n"
            "##OctopuSV_mode=caller\n"
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n"
            "chr1\t100\tbad\tN\t<DEL>\t60\tPASS\t"
            "SVTYPE=DEL;END=oops;SVLEN=50;CHR2=chr1;SUPPORT=.;SVMETHOD=OctopuSV;"
            "RTID=.;AF=.;STRAND=.;RNAMES=.\t"
            "GT:AD:LN:ST:QV:TY:ID:SC:REF:ALT:CO\t"
            "0/1:.,.:50:.:60:DEL:bad:TestCaller:N:<DEL>:chr1_100-chr1_oops\n"
        )

    monkeypatch.setattr("octopusv.cli.convert.write_sv_vcf", invalid_but_successful_writer)

    result = _invoke_correct(input_vcf, output_svcf)

    assert result.exit_code != 0
    assert "validation" in result.output.lower()
    assert output_svcf.read_text() == "ORIGINAL\n"
    assert not list(tmp_path.glob(".out.svcf.*.tmp"))


def test_correct_does_not_treat_format_dr_as_alt_support(tmp_path):
    input_vcf = tmp_path / "dr_only.vcf"
    output_svcf = tmp_path / "out.svcf"
    input_vcf.write_text(
        _header(
            fmt_lines=[
                '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
                '##FORMAT=<ID=DR,Number=1,Type=Integer,Description="Reference reads">',
            ]
        )
        + _record(
            record_id="dr_only",
            info="SVTYPE=DEL;END=150;SVLEN=-50",
            fmt="GT:DR",
            sample="0/1:50",
        )
    )

    # DR is reference-read depth, not ALT support.  A max-support filter must
    # therefore not discard this event merely because DR is large.
    result = _invoke_correct(input_vcf, output_svcf, "--max-support", "10")

    assert result.exit_code == 0, result.output
    records = _data_lines(output_svcf)
    assert len(records) == 1
    assert records[0].split("\t")[2] == "dr_only"


def test_correct_allows_truly_empty_vcf_and_writes_valid_empty_svcf(tmp_path):
    input_vcf = tmp_path / "empty.vcf"
    output_svcf = tmp_path / "empty.svcf"
    input_vcf.write_text(_header())

    result = _invoke_correct(input_vcf, output_svcf)

    assert result.exit_code == 0, result.output
    assert output_svcf.exists()
    assert _data_lines(output_svcf) == []
    assert "no data records" in result.output.lower()
    _assert_valid_svcf(output_svcf)
