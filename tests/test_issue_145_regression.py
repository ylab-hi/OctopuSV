from __future__ import annotations

from pathlib import Path

import pytest
from typer.testing import CliRunner

pytest.importorskip("natsort")

from octopusv.cli.cli import app
from octopusv.utils.svcf_validator import SVCFValidator


runner = CliRunner()


def _header() -> str:
    return (
        "##fileformat=VCFv4.2\n"
        "##source=DRAGEN\n"
        "##contig=<ID=chr1,length=1000000>\n"
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n"
    )


def _cnv_record(
    *,
    pos: int,
    record_id: str,
    alt: str,
    end: int,
    svlen: int,
) -> str:
    return (
        f"chr1\t{pos}\t{record_id}\tN\t{alt}\t60\tPASS\t"
        f"SVTYPE=CNV;END={end};SVLEN={svlen}\tGT\t0/1\n"
    )


def _info_dict(record_line: str) -> dict[str, str]:
    info = record_line.split("\t")[7]
    return {
        key: value
        for item in info.split(";")
        if "=" in item
        for key, value in [item.split("=", 1)]
    }


def test_issue_145_dragen_cnv_is_normalized_to_supported_svtypes(tmp_path: Path):
    """DRAGEN CNV records must become standard DEL/DUP SVCF events."""
    input_vcf = tmp_path / "dragen_cnv.vcf"
    output_svcf = tmp_path / "dragen_cnv.svcf"

    input_vcf.write_text(
        _header()
        + _cnv_record(
            pos=100,
            record_id="dragen_loss",
            alt="<DEL>",
            end=150,
            svlen=-50,
        )
        + _cnv_record(
            pos=300,
            record_id="dragen_gain",
            alt="<DUP>",
            end=400,
            svlen=100,
        ),
        encoding="utf-8",
    )

    result = runner.invoke(
        app,
        [
            "correct",
            "-i",
            str(input_vcf),
            "-o",
            str(output_svcf),
        ],
    )

    assert result.exit_code == 0, (
        f"correct failed:\n{result.output}\nexception: {result.exception!r}"
    )
    assert output_svcf.exists()

    records = [
        line
        for line in output_svcf.read_text(encoding="utf-8").splitlines()
        if line and not line.startswith("#")
    ]
    assert len(records) == 2

    by_id = {line.split("\t")[2]: _info_dict(line) for line in records}

    assert by_id["dragen_loss"]["SVTYPE"] == "DEL"
    assert by_id["dragen_gain"]["SVTYPE"] == "DUP"

    # The unsupported intermediate CNV type must not leak into SVCF 1.1.
    assert all(info["SVTYPE"] != "CNV" for info in by_id.values())

    validator = SVCFValidator(str(output_svcf))
    validator.validate()
    assert validator.errors == [], validator.to_summary(max_issues=None)
