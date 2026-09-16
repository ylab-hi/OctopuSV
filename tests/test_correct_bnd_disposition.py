from __future__ import annotations

from pathlib import Path

import pytest
from typer.testing import CliRunner

pytest.importorskip("natsort")

from octopusv.cli.cli import app


runner = CliRunner()


def _bnd(chrom: str, pos: int, record_id: str, alt: str) -> str:
    return (
        f"{chrom}\t{pos}\t{record_id}\tN\t{alt}\t60\tPASS\t"
        "SVTYPE=BND\tGT\t0/1\n"
    )


def _del(chrom: str, pos: int, record_id: str, end: int) -> str:
    svlen = end - pos
    return (
        f"{chrom}\t{pos}\t{record_id}\tN\t<DEL>\t60\tPASS\t"
        f"SVTYPE=DEL;END={end};SVLEN=-{svlen}\tGT\t0/1\n"
    )


def _raw_vcf(records: list[str], *, extra_contigs: list[str] | None = None) -> str:
    contigs = ["chr1", "chr2"] + list(extra_contigs or [])
    return (
        "##fileformat=VCFv4.2\n"
        "##source=DispositionRegressionCaller\n"
        + "".join(f"##contig=<ID={contig},length=1000000>\n" for contig in contigs)
        + "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n"
        + "".join(records)
    )


def _invoke(tmp_path: Path, records: list[str], *extra_args: str, extra_contigs=None):
    input_vcf = tmp_path / "input.vcf"
    output_svcf = tmp_path / "output.svcf"
    input_vcf.write_text(_raw_vcf(records, extra_contigs=extra_contigs))
    result = runner.invoke(
        app,
        ["correct", "-i", str(input_vcf), "-o", str(output_svcf), *extra_args],
    )
    return result, output_svcf


def _records(path: Path) -> list[list[str]]:
    if not path.exists():
        return []
    return [
        line.split("\t")
        for line in path.read_text().splitlines()
        if line and not line.startswith("#")
    ]


def _info(fields: list[str]) -> dict[str, str]:
    return dict(item.split("=", 1) for item in fields[7].split(";"))


def test_r2_unhandled_same_direction_special_pair_falls_back_to_single_tra(tmp_path):
    result, output = _invoke(
        tmp_path,
        [
            _bnd("chr1", 100, "r2.a", "N[chr2:200["),
            _bnd("chr1", 102, "r2.b", "N[chr2:201["),
        ],
    )

    assert result.exit_code == 0, result.output
    rows = _records(output)
    assert [row[2] for row in rows] == ["r2.a", "r2.b"]
    for row in rows:
        info = _info(row)
        assert info["SVTYPE"] == "TRA"
        assert info["CHR2"] == "chr2"
        assert info["SVLEN"] == "."


def test_r3_duplicate_coordinate_record_not_lost_when_other_copy_is_mate_paired(tmp_path):
    result, output = _invoke(
        tmp_path,
        [
            _bnd("chr1", 100, "r3.a", "N]chr2:200]"),
            _bnd("chr1", 100, "r3.aprime", "N]chr2:200]"),
            _bnd("chr2", 200, "r3.b", "[chr1:100[N"),
        ],
    )

    assert result.exit_code == 0, result.output
    rows = _records(output)
    ids = [row[2] for row in rows]
    assert sorted(ids) == ["r3.a", "r3.aprime", "r3.b"]
    assert len(ids) == len(set(ids)) == 3
    assert all(_info(row)["SVTYPE"] == "TRA" for row in rows)


def test_r4_unrecognized_iupac_bnd_alt_fails_instead_of_disappearing(tmp_path):
    result, output = _invoke(
        tmp_path,
        [_bnd("chr1", 100, "r4.bad", "R[chr2:200[")],
    )

    assert result.exit_code != 0
    assert "Malformed or unsupported BND ALT" in result.output
    assert "r4.bad" in result.output
    assert not output.exists()


def test_r5_true_single_breakend_fails_by_default_with_actionable_summary(tmp_path):
    result, output = _invoke(
        tmp_path,
        [
            _bnd("chr1", 100, "r5.a", "N."),
            _bnd("chr1", 200, "r5.b", ".N"),
        ],
    )

    assert result.exit_code != 0
    assert "Found 2 true single-breakend BND record(s)" in result.output
    assert "--skip-single-breakends" in result.output
    assert "r5.a" in result.output
    assert "r5.b" in result.output
    assert not output.exists()


def test_r5_explicit_skip_is_counted_in_output_header_and_other_records_survive(tmp_path):
    result, output = _invoke(
        tmp_path,
        [
            _bnd("chr1", 100, "r5.single", "N."),
            _del("chr1", 300, "kept.del", 350),
        ],
        "--skip-single-breakends",
    )

    assert result.exit_code == 0, result.output
    text = output.read_text()
    assert "##OctopuSV_skipped_single_breakends=1" in text
    rows = _records(output)
    assert [row[2] for row in rows] == ["kept.del"]
    assert _info(rows[0])["SVTYPE"] == "DEL"
    assert "Skipped 1 true single-breakend record(s)" in result.output

    validation = runner.invoke(app, ["validate-svcf", str(output)])
    assert validation.exit_code == 0, validation.output
    assert "Status: PASS" in validation.output


def test_skip_single_breakends_does_not_hide_malformed_bnd_alt(tmp_path):
    result, output = _invoke(
        tmp_path,
        [_bnd("chr1", 100, "still.bad", "R[chr2:200[")],
        "--skip-single-breakends",
    )

    assert result.exit_code != 0
    assert "Malformed or unsupported BND ALT" in result.output
    assert not output.exists()


def test_r6_colon_contig_is_explicitly_rejected_until_full_stack_support_exists(tmp_path):
    hla = "HLA-A*01:01:01:01"
    result, output = _invoke(
        tmp_path,
        [_bnd("chr1", 100, "r6.hla", f"N[{hla}:200[")],
        extra_contigs=[hla],
    )

    assert result.exit_code != 0
    assert "Unsupported BND remote contig" in result.output
    assert "contig names containing ':'" in result.output
    assert "r6.hla" in result.output
    assert not output.exists()


@pytest.mark.parametrize(
    ("record", "extra_contigs"),
    [
        (
            _del("HLA-A*01:01:01:01", 100, "colon.local", 150),
            ["HLA-A*01:01:01:01"],
        ),
        (
            (
                "chr1\t100\tcolon.chr2\tN\t<DEL>\t60\tPASS\t"
                "SVTYPE=DEL;END=150;SVLEN=-50;CHR2=HLA-A*01:01:01:01\tGT\t0/1\n"
            ),
            ["HLA-A*01:01:01:01"],
        ),
    ],
)
def test_colon_contig_is_rejected_before_svcf_serialization(tmp_path, record, extra_contigs):
    result, output = _invoke(
        tmp_path,
        [record],
        extra_contigs=extra_contigs,
    )

    assert result.exit_code != 0
    assert "Unsupported contig name" in result.output
    assert "contig names containing ':'" in result.output
    assert not output.exists()
