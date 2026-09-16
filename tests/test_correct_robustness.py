from __future__ import annotations

from pathlib import Path

import pytest
from typer.testing import CliRunner

pytest.importorskip("natsort")

from octopusv.cli.cli import app
from octopusv.utils.normal_vcf_parser import parse_vcf


runner = CliRunner()


def _invoke_correct(input_vcf: Path, output_svcf: Path):
    return runner.invoke(
        app,
        ["correct", "-i", str(input_vcf), "-o", str(output_svcf)],
    )


def _header(sample_names: list[str], *, source: str = "TestCaller") -> str:
    return (
        "##fileformat=VCFv4.2\n"
        f"##source={source}\n"
        "##contig=<ID=chr1,length=1000000>\n"
        "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n"
        + "\t".join(
            [
                "#CHROM",
                "POS",
                "ID",
                "REF",
                "ALT",
                "QUAL",
                "FILTER",
                "INFO",
                "FORMAT",
                *sample_names,
            ]
        )
        + "\n"
    )


def _del_record(sample_values: list[str]) -> str:
    return (
        "chr1\t100\tr1\tN\t<DEL>\t60\tPASS\t"
        "SVTYPE=DEL;END=150;SVLEN=-50\tGT\t"
        + "\t".join(sample_values)
        + "\n"
    )


def test_correct_accepts_ordinary_four_sample_vcf(tmp_path):
    input_vcf = tmp_path / "four_samples.vcf"
    output_svcf = tmp_path / "four_samples.svcf"
    sample_names = ["S1", "S2", "S3", "S4"]
    sample_values = ["0/1", "0/0", "1/1", "./."]
    input_vcf.write_text(_header(sample_names) + _del_record(sample_values))

    result = _invoke_correct(input_vcf, output_svcf)

    assert result.exit_code == 0, result.output
    lines = output_svcf.read_text().splitlines()
    chrom_header = next(line for line in lines if line.startswith("#CHROM"))
    record = next(line for line in lines if line and not line.startswith("#"))
    assert chrom_header.split("\t")[9:] == sample_names
    assert len(record.split("\t")) == 9 + len(sample_names)


def test_svaba_13_column_special_path_is_preserved(tmp_path):
    input_vcf = tmp_path / "svaba.vcf"
    output_svcf = tmp_path / "svaba.svcf"
    sample_names = ["S1", "S2", "S3", "S4"]
    sample_values = ["0/0", "0/1", "1/1", "./."]
    input_vcf.write_text(
        _header(sample_names, source="SvABA") + _del_record(sample_values)
    )

    _, _, _, non_bnd = parse_vcf(input_vcf)
    assert len(non_bnd) == 1
    # Historical SvABA behavior is intentionally preserved: the relevant
    # sample is the final (13th) column, not a generic four-sample expansion.
    assert non_bnd[0].samples == [sample_values[-1]]
    assert non_bnd[0].sample == sample_values[-1]

    result = _invoke_correct(input_vcf, output_svcf)
    assert result.exit_code == 0, result.output
    lines = output_svcf.read_text().splitlines()
    chrom_header = next(line for line in lines if line.startswith("#CHROM"))
    record = next(line for line in lines if line and not line.startswith("#"))
    assert chrom_header.split("\t")[9:] == [sample_names[-1]]
    assert len(record.split("\t")) == 10


def test_correct_rejects_record_width_mismatch(tmp_path, caplog):
    input_vcf = tmp_path / "truncated.vcf"
    output_svcf = tmp_path / "out.svcf"
    input_vcf.write_text(
        _header(["S1", "S2"]) + _del_record(["0/1"])
    )

    result = _invoke_correct(input_vcf, output_svcf)

    assert result.exit_code != 0
    assert "#CHROM declares 11 columns, but the record has 10" in caplog.text
    assert not output_svcf.exists()


def test_correct_requires_real_chrom_header(tmp_path, caplog):
    input_vcf = tmp_path / "missing_chrom.vcf"
    output_svcf = tmp_path / "out.svcf"
    input_vcf.write_text(
        "##fileformat=VCFv4.2\n"
        "##source=TestCaller\n"
        "##contig=<ID=chr1,length=1000000>\n"
        "chr1\t100\tr1\tN\t<DEL>\t60\tPASS\t"
        "SVTYPE=DEL;END=150;SVLEN=-50\tGT\t0/1\n"
    )

    result = _invoke_correct(input_vcf, output_svcf)

    assert result.exit_code != 0
    assert "before a #CHROM header" in caplog.text or "must contain a #CHROM" in caplog.text
    assert not output_svcf.exists()


def test_correct_writer_failure_preserves_existing_output(tmp_path, monkeypatch):
    input_vcf = tmp_path / "input.vcf"
    output_svcf = tmp_path / "out.svcf"
    input_vcf.write_text(_header(["S1"]) + _del_record(["0/1"]))
    output_svcf.write_text("ORIGINAL\n")

    def failing_writer(contig_lines, events, output_file, input_vcf_file=None):
        Path(output_file).write_text("PARTIAL\n")
        raise RuntimeError("simulated writer failure")

    monkeypatch.setattr("octopusv.cli.convert.write_sv_vcf", failing_writer)

    result = _invoke_correct(input_vcf, output_svcf)

    assert result.exit_code != 0
    assert output_svcf.read_text() == "ORIGINAL\n"
    assert not list(tmp_path.glob(".out.svcf.*.tmp"))


def test_correct_atomic_write_replaces_existing_output_on_success(tmp_path):
    input_vcf = tmp_path / "input.vcf"
    output_svcf = tmp_path / "out.svcf"
    input_vcf.write_text(_header(["S1"]) + _del_record(["0/1"]))
    output_svcf.write_text("OLD\n")

    result = _invoke_correct(input_vcf, output_svcf)

    assert result.exit_code == 0, result.output
    assert output_svcf.read_text() != "OLD\n"
    assert "##SVCFVersion=1.1" in output_svcf.read_text()
    assert not list(tmp_path.glob(".out.svcf.*.tmp"))
