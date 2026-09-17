from __future__ import annotations

from pathlib import Path

import pytest
from typer.testing import CliRunner

pytest.importorskip("natsort")

from octopusv.cli.cli import app
from octopusv.utils.sample_consensus import resolve_sample_consensus
from octopusv.utils.svcf_sample_parser import parse_svcf_sample_block


RUNNER = CliRunner()


def _record_line(*, fmt: str | None = None, sample: str | None = None, record_id: str = "sv1", pos: int = 100) -> str:
    end = pos + 50
    fields = [
        "chr1",
        str(pos),
        record_id,
        "N",
        "<DEL>",
        "60",
        "PASS",
        f"SVTYPE=DEL;END={end};SVLEN=-50",
    ]
    if fmt is not None:
        fields.append(fmt)
    if sample is not None:
        fields.append(sample)
    return "\t".join(fields) + "\n"


def _write_vcf(path: Path, records: list[str], *, fmt_header: str = "") -> None:
    # fmt_header is the literal suffix appended after INFO in the #CHROM line,
    # e.g. "\tFORMAT" for a 9-column site-only VCF or "\tFORMAT\tS1" for
    # an ordinary single-sample VCF.
    path.write_text(
        "##fileformat=VCFv4.2\n"
        "##source=Step3Test\n"
        "##contig=<ID=chr1,length=1000000>\n"
        f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO{fmt_header}\n"
        + "".join(records)
    )


def _run_correct(input_vcf: Path, output_svcf: Path):
    return RUNNER.invoke(
        app,
        ["correct", "-i", str(input_vcf), "-o", str(output_svcf)],
    )


def _first_record(path: Path) -> list[str]:
    with path.open() as handle:
        for line in handle:
            if not line.startswith("#"):
                return line.rstrip("\n").split("\t")
    raise AssertionError("No output record found")


def _records(path: Path) -> list[list[str]]:
    result = []
    with path.open() as handle:
        for line in handle:
            if not line.startswith("#"):
                result.append(line.rstrip("\n").split("\t"))
    return result


def _sample(record: list[str]) -> dict[str, str]:
    return parse_svcf_sample_block(record[8], record[9])


def test_eight_column_site_only_vcf_becomes_carrier_unknown(tmp_path):
    raw = tmp_path / "site_only_8.vcf"
    out = tmp_path / "site_only_8.svcf"
    _write_vcf(raw, [_record_line()])

    result = _run_correct(raw, out)

    assert result.exit_code == 0, result.output
    sample = _sample(_first_record(out))
    assert sample["GT"] == "1/."
    assert sample["AD"] == ".,."
    assert "all explicit source GT values are missing" not in result.output


@pytest.mark.parametrize(
    ("raw_format", "expected_gt", "expected_ad"),
    [
        ("GT:DP", "1/.", ".,."),
        ("DP:AD", "1/.", ".,."),
        (".", "1/.", ".,."),
    ],
)
def test_nine_column_site_only_vcf_synthesizes_aligned_carrier_unknown_sample(
    tmp_path,
    raw_format,
    expected_gt,
    expected_ad,
):
    raw = tmp_path / "site_only_9.vcf"
    out = tmp_path / "site_only_9.svcf"
    _write_vcf(raw, [_record_line(fmt=raw_format)], fmt_header="\tFORMAT")

    result = _run_correct(raw, out)

    assert result.exit_code == 0, result.output
    sample = _sample(_first_record(out))
    assert sample["GT"] == expected_gt
    assert sample["AD"] == expected_ad
    assert "all explicit source GT values are missing" not in result.output


def test_explicit_missing_gt_is_preserved_and_warned_once(tmp_path):
    raw = tmp_path / "explicit_missing.vcf"
    out = tmp_path / "explicit_missing.svcf"
    _write_vcf(
        raw,
        [
            _record_line(fmt="GT:DP", sample="./.:20", record_id="sv1", pos=100),
            _record_line(fmt="GT:DP", sample=".:30", record_id="sv2", pos=300),
        ],
        fmt_header="\tFORMAT\tS1",
    )

    result = _run_correct(raw, out)

    assert result.exit_code == 0, result.output
    records = _records(out)
    assert [_sample(record)["GT"] for record in records] == ["./.", "."]
    assert result.output.count("all explicit source GT values are missing") == 1
    assert "will not contribute carrier/absence votes" in result.output


def test_explicit_informative_gt_suppresses_all_missing_warning(tmp_path):
    raw = tmp_path / "mixed_gt.vcf"
    out = tmp_path / "mixed_gt.svcf"
    _write_vcf(
        raw,
        [
            _record_line(fmt="GT", sample="./.", record_id="sv1", pos=100),
            _record_line(fmt="GT", sample="0/1", record_id="sv2", pos=300),
        ],
        fmt_header="\tFORMAT\tS1",
    )

    result = _run_correct(raw, out)

    assert result.exit_code == 0, result.output
    assert [_sample(record)["GT"] for record in _records(out)] == ["./.", "0/1"]
    assert "all explicit source GT values are missing" not in result.output


def test_site_only_carrier_unknown_and_explicit_het_preserve_presence_but_not_zygosity():
    consensus = resolve_sample_consensus(
        [
            ("site_only", "1/."),
            ("genotyped", "0/1"),
        ]
    )

    assert consensus.gt == "1/."
    assert consensus.unique_carrier_callers == 2
    assert consensus.valid_caller_votes == 2


def test_explicit_missing_gt_remains_no_vote_in_consensus():
    consensus = resolve_sample_consensus([("caller", "./.")])

    assert consensus.gt == "./."
    assert consensus.unique_carrier_callers == 0
    assert consensus.valid_caller_votes == 0
