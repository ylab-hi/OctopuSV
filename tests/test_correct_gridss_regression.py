from __future__ import annotations

from pathlib import Path

import pytest
from typer.testing import CliRunner

pytest.importorskip("natsort")

from octopusv.cli.cli import app
from octopusv.utils.svcf_sample_parser import parse_svcf_sample_block
from octopusv.utils.svcf_validator import SVCFValidator


RUNNER = CliRunner()


def _gridss_header(*, contigs: tuple[str, ...] = ("chr1", "chr2")) -> str:
    """Minimal GRIDSS-like single-sample VCF header used by these regressions."""
    return (
        "##fileformat=VCFv4.2\n"
        "##source=GRIDSS\n"
        + "".join(
            f"##contig=<ID={contig},length=1000000>\n" for contig in contigs
        )
        + '##FILTER=<ID=LOW_QUAL,Description="Low quality GRIDSS call">\n'
        + '##INFO=<ID=MATEID,Number=1,Type=String,Description="Mate breakend ID">\n'
        + '##INFO=<ID=EVENT,Number=1,Type=String,Description="Breakend event ID">\n'
        + '##INFO=<ID=CIPOS,Number=2,Type=Integer,Description="Confidence interval around POS">\n'
        + '##INFO=<ID=HOMLEN,Number=1,Type=Integer,Description="Microhomology length">\n'
        + '##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description="Imprecise breakpoint">\n'
        + '##INFO=<ID=BEALN,Number=1,Type=String,Description="Breakend alignment annotation">\n'
        + '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
        + "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n"
    )


def _bnd(
    chrom: str,
    pos: int,
    record_id: str,
    alt: str,
    *,
    qual: float = 60.0,
    filt: str = "PASS",
    extra_info: str = "",
    gt: str = ".",
) -> str:
    info = "SVTYPE=BND"
    if extra_info:
        info += ";" + extra_info
    return (
        f"{chrom}\t{pos}\t{record_id}\tN\t{alt}\t{qual}\t{filt}\t"
        f"{info}\tGT\t{gt}\n"
    )


def _run_correct(
    tmp_path: Path,
    records: list[str],
    *extra_args: str,
    contigs: tuple[str, ...] = ("chr1", "chr2"),
):
    input_vcf = tmp_path / "gridss.vcf"
    output_svcf = tmp_path / "gridss.svcf"
    input_vcf.write_text(_gridss_header(contigs=contigs) + "".join(records))
    result = RUNNER.invoke(
        app,
        ["correct", "-i", str(input_vcf), "-o", str(output_svcf), *extra_args],
    )
    return result, output_svcf


def _records(path: Path) -> list[list[str]]:
    if not path.exists():
        return []
    return [
        line.rstrip("\n").split("\t")
        for line in path.read_text().splitlines()
        if line and not line.startswith("#")
    ]


def _info(row: list[str]) -> dict[str, str]:
    return {
        key: value
        for item in row[7].split(";")
        if "=" in item
        for key, value in [item.split("=", 1)]
    }


def _sample(row: list[str]) -> dict[str, str]:
    return parse_svcf_sample_block(row[8], row[9])


def _assert_valid_svcf(path: Path) -> None:
    validator = SVCFValidator(str(path))
    validator.validate()
    assert validator.errors == [], validator.to_summary(max_issues=None)


def test_gridss_interchromosomal_mate_pair_with_inserted_sequence_becomes_tra(tmp_path):
    """GRIDSS-style bracket ALT may carry inserted sequence before the mate target."""
    result, output = _run_correct(
        tmp_path,
        [
            _bnd(
                "chr1",
                1000,
                "gridss.tra.1",
                "AGGT[chr2:5000[",
                qual=123.4,
                extra_info=(
                    "MATEID=gridss.tra.2;EVENT=event.tra;CIPOS=-2,2;"
                    "HOMLEN=0;IMPRECISE;BEALN=chr1:995-1004|chr2:4995-5004"
                ),
            ),
            _bnd(
                "chr2",
                5000,
                "gridss.tra.2",
                "]chr1:1000]T",
                qual=99.5,
                extra_info="MATEID=gridss.tra.1;EVENT=event.tra",
            ),
        ],
    )

    assert result.exit_code == 0, result.output
    rows = _records(output)
    assert len(rows) == 1

    row = rows[0]
    info = _info(row)
    assert info["SVTYPE"] == "TRA"
    assert row[0] == "chr1"
    assert row[1] == "1000"
    assert info["CHR2"] == "chr2"
    assert info["END"] == "5000"
    assert info["SVLEN"] == "."

    # GRIDSS explicitly reports missing genotype. Preserve that source fact and
    # source identity instead of synthesizing carrier/zygosity information.
    sample = _sample(row)
    assert sample["GT"] == "."
    assert sample["SC"] == "GRIDSS"
    assert sample["ALT"] == "AGGT[chr2:5000["
    assert "all explicit source GT values are missing" in result.output
    _assert_valid_svcf(output)


def test_gridss_same_chromosome_del_geometry_becomes_del(tmp_path):
    result, output = _run_correct(
        tmp_path,
        [
            _bnd(
                "chr1",
                10000,
                "gridss.del.1",
                "AG[chr1:30000[",
                qual=50.5,
                extra_info="MATEID=gridss.del.2;EVENT=event.del",
            ),
            _bnd(
                "chr1",
                30000,
                "gridss.del.2",
                "]chr1:10000]T",
                qual=49.0,
                extra_info="MATEID=gridss.del.1;EVENT=event.del",
            ),
        ],
        contigs=("chr1",),
    )

    assert result.exit_code == 0, result.output
    rows = _records(output)
    assert len(rows) == 1
    row = rows[0]
    info = _info(row)
    assert row[4] == "<DEL>"
    assert info["SVTYPE"] == "DEL"
    assert row[1] == "10000"
    assert info["END"] == "30000"
    assert info["SVLEN"] == "20000"
    assert info["CHR2"] == "chr1"
    _assert_valid_svcf(output)


def test_gridss_same_chromosome_inv_geometry_preserves_nonpass_filter(tmp_path):
    result, output = _run_correct(
        tmp_path,
        [
            _bnd(
                "chr1",
                40000,
                "gridss.inv.1",
                "AC]chr1:60000]",
                qual=31.2,
                filt="LOW_QUAL",
                extra_info="MATEID=gridss.inv.2;EVENT=event.inv;IMPRECISE",
            ),
            _bnd(
                "chr1",
                60000,
                "gridss.inv.2",
                "GT]chr1:40000]",
                qual=30.1,
                filt="LOW_QUAL",
                extra_info="MATEID=gridss.inv.1;EVENT=event.inv;IMPRECISE",
            ),
        ],
        contigs=("chr1",),
    )

    assert result.exit_code == 0, result.output
    rows = _records(output)
    assert len(rows) == 1
    row = rows[0]
    info = _info(row)
    assert row[4] == "<INV>"
    assert row[6] == "LOW_QUAL"
    assert info["SVTYPE"] == "INV"
    assert info["END"] == "60000"
    assert info["SVLEN"] == "20000"
    _assert_valid_svcf(output)


def test_gridss_true_single_breakend_fails_by_default_with_actionable_message(tmp_path):
    result, output = _run_correct(
        tmp_path,
        [
            _bnd(
                "chr1",
                70000,
                "gridss.single.1",
                "GAATTC.",
                qual=42.0,
                extra_info="EVENT=event.single",
            )
        ],
        contigs=("chr1",),
    )

    assert result.exit_code != 0
    assert "true single-breakend BND record(s)" in result.output
    assert "gridss.single.1" in result.output
    assert "--skip-single-breakends" in result.output
    assert not output.exists()


def test_gridss_single_breakend_can_be_explicitly_skipped_and_output_validates(tmp_path):
    result, output = _run_correct(
        tmp_path,
        [
            _bnd(
                "chr1",
                70000,
                "gridss.single.1",
                "GAATTC.",
                qual=42.0,
                extra_info="EVENT=event.single",
            )
        ],
        "--skip-single-breakends",
        contigs=("chr1",),
    )

    assert result.exit_code == 0, result.output
    text = output.read_text()
    assert "##OctopuSV_skipped_single_breakends=1" in text
    assert _records(output) == []
    assert "Skipped 1 true single-breakend record(s)" in result.output

    validation = RUNNER.invoke(app, ["validate-svcf", str(output)])
    assert validation.exit_code == 0, validation.output
    assert "Status: PASS" in validation.output


def test_gridss_explicit_missing_gt_is_not_rewritten_to_carrier_unknown(tmp_path):
    result, output = _run_correct(
        tmp_path,
        [
            _bnd(
                "chr1",
                80000,
                "gridss.gt.missing",
                "ACGT[chr2:90000[",
                qual=80.0,
                extra_info="EVENT=event.gt",
                gt=".",
            )
        ],
    )

    assert result.exit_code == 0, result.output
    row = _records(output)[0]
    sample = _sample(row)
    assert sample["GT"] == "."
    assert sample["GT"] != "1/."
    assert "all explicit source GT values are missing" in result.output
    _assert_valid_svcf(output)
