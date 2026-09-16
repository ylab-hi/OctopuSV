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


def _sv(
    chrom: str,
    pos: int,
    record_id: str,
    alt: str,
    svtype: str,
    *,
    end: int | str,
    svlen: int | str,
) -> str:
    return (
        f"{chrom}\t{pos}\t{record_id}\tN\t{alt}\t60\tPASS\t"
        f"SVTYPE={svtype};END={end};SVLEN={svlen}\tGT\t0/1\n"
    )


def _raw_vcf(records: list[str]) -> str:
    return (
        "##fileformat=VCFv4.2\n"
        "##source=RegressionCaller\n"
        "##contig=<ID=chr1,length=1000000>\n"
        "##contig=<ID=chr2,length=1000000>\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n"
        + "".join(records)
    )


def _run_correct(tmp_path: Path, records: list[str]) -> list[dict[str, str]]:
    input_vcf = tmp_path / "input.vcf"
    output_svcf = tmp_path / "output.svcf"
    input_vcf.write_text(_raw_vcf(records))

    result = runner.invoke(
        app,
        ["correct", "-i", str(input_vcf), "-o", str(output_svcf)],
    )
    assert result.exit_code == 0, result.output

    parsed: list[dict[str, str]] = []
    for line in output_svcf.read_text().splitlines():
        if not line or line.startswith("#"):
            continue
        fields = line.split("\t")
        info = {}
        for item in fields[7].split(";"):
            key, value = item.split("=", 1)
            info[key] = value
        parsed.append(
            {
                "CHROM": fields[0],
                "POS": fields[1],
                "ID": fields[2],
                "ALT": fields[4],
                **info,
            }
        )
    return parsed


def _single(records: list[dict[str, str]]) -> dict[str, str]:
    assert len(records) == 1
    return records[0]


def test_correct_same_chr_bnd_pair_to_del_is_locked(tmp_path):
    records = _run_correct(
        tmp_path,
        [
            _bnd("chr1", 100, "del.left", "N[chr1:200["),
            _bnd("chr1", 200, "del.right", "]chr1:100]N"),
        ],
    )

    event = _single(records)
    assert event["SVTYPE"] == "DEL"
    assert event["POS"] == "100"
    assert event["END"] == "200"
    assert event["SVLEN"] == "100"
    assert event["CHR2"] == "chr1"
    assert event["ALT"] == "<DEL>"


def test_correct_del_geometry_with_long_inserted_sequence_to_ins_is_locked(tmp_path):
    inserted = "A" * 80
    records = _run_correct(
        tmp_path,
        [
            _bnd("chr1", 100, "ins.left", f"N{inserted}[chr1:200["),
            _bnd("chr1", 200, "ins.right", "]chr1:100]N"),
        ],
    )

    event = _single(records)
    assert event["SVTYPE"] == "INS"
    assert event["POS"] == "100"
    assert event["SVLEN"] == "80"
    # Locked SVCF 1.1 internal insertion convention.
    assert event["END"] == "180"
    assert event["CHR2"] == "chr1"
    assert event["ALT"] == "<INS>"


def test_correct_same_chr_bnd_pair_to_dup_is_locked(tmp_path):
    records = _run_correct(
        tmp_path,
        [
            _bnd("chr1", 200, "dup.right", "N[chr1:100["),
            _bnd("chr1", 100, "dup.left", "]chr1:200]N"),
        ],
    )

    event = _single(records)
    assert event["SVTYPE"] == "DUP"
    assert event["POS"] == "100"
    assert event["END"] == "200"
    assert event["SVLEN"] == "100"
    assert event["CHR2"] == "chr1"
    assert event["ALT"] == "<DUP>"


def test_correct_same_chr_bnd_pair_to_inv_is_locked(tmp_path):
    records = _run_correct(
        tmp_path,
        [
            _bnd("chr1", 100, "inv.left", "N]chr1:200]"),
            _bnd("chr1", 200, "inv.right", "N]chr1:100]"),
        ],
    )

    event = _single(records)
    assert event["SVTYPE"] == "INV"
    assert event["POS"] == "100"
    assert event["END"] == "200"
    assert event["SVLEN"] == "100"
    assert event["CHR2"] == "chr1"
    assert event["ALT"] == "<INV>"


def test_correct_unresolved_same_chr_bnd_with_explicit_mate_stays_bnd(tmp_path):
    event = _single(
        _run_correct(
            tmp_path,
            [_bnd("chr1", 100, "kept.bnd", "N[chr1:300[")],
        )
    )

    assert event["SVTYPE"] == "BND"
    assert event["POS"] == "100"
    assert event["END"] == "300"
    assert event["CHR2"] == "chr1"
    assert event["SVLEN"] == "."
    assert event["ALT"] == "N[chr1:300["


def test_correct_single_interchrom_bnd_with_explicit_mate_becomes_tra(tmp_path):
    """Intentional OctopuSV design: a remote mate coordinate is sufficient for TRA."""
    event = _single(
        _run_correct(
            tmp_path,
            [_bnd("chr1", 100, "single.tra", "N[chr2:200[")],
        )
    )

    assert event["SVTYPE"] == "TRA"
    assert event["POS"] == "100"
    assert event["END"] == "200"
    assert event["CHR2"] == "chr2"
    assert event["SVLEN"] == "."
    assert event["ALT"] == "N[chr2:200["


def test_correct_reciprocal_interchrom_pair_to_tra_is_locked(tmp_path):
    records = _run_correct(
        tmp_path,
        [
            _bnd("chr1", 100, "tra.a", "N]chr2:200]"),
            _bnd("chr2", 200, "tra.b", "[chr1:100[N"),
        ],
    )

    assert len(records) == 2
    by_id = {record["ID"]: record for record in records}
    assert set(by_id) == {"tra.a", "tra.b"}
    assert by_id["tra.a"]["SVTYPE"] == "TRA"
    assert by_id["tra.a"]["CHR2"] == "chr2"
    assert by_id["tra.a"]["END"] == "200"
    assert by_id["tra.a"]["RTID"] == "tra.b"
    assert by_id["tra.b"]["SVTYPE"] == "TRA"
    assert by_id["tra.b"]["CHR2"] == "chr1"
    assert by_id["tra.b"]["END"] == "100"
    assert by_id["tra.b"]["RTID"] == "tra.a"


def test_correct_direction_different_special_pair_to_tra_is_locked(tmp_path):
    records = _run_correct(
        tmp_path,
        [
            _bnd("chr1", 100, "special.a", "N[chr2:200["),
            _bnd("chr1", 102, "special.b", "N]chr2:201]"),
        ],
    )

    assert len(records) == 2
    by_id = {record["ID"]: record for record in records}
    assert set(by_id) == {"special.a", "special.b"}
    for record in by_id.values():
        assert record["SVTYPE"] == "TRA"
        assert record["CHR2"] == "chr2"
        assert record["SVLEN"] == "."


def test_correct_native_ins_uses_locked_svcf_internal_endpoint(tmp_path):
    event = _single(
        _run_correct(
            tmp_path,
            [
                _sv(
                    "chr1",
                    100,
                    "native.ins",
                    "<INS>",
                    "INS",
                    end=100,
                    svlen=47,
                )
            ],
        )
    )

    assert event["SVTYPE"] == "INS"
    assert event["POS"] == "100"
    assert event["SVLEN"] == "47"
    assert event["END"] == "147"
    assert event["CHR2"] == "chr1"


def test_characterize_pending_decision_overlapping_mate_pair_strategies(tmp_path):
    """Characterization only: do not treat this current duplicate as desired science.

    The {t[p[, t]p]} mate pair currently matches both the independent and merge
    TRA converters.  Keep this test visibly named as a pending-decision snapshot
    until R1 is resolved with stronger semantic evidence.
    """
    records = _run_correct(
        tmp_path,
        [
            _bnd("chr1", 100, "overlap.a", "N[chr2:200["),
            _bnd("chr2", 200, "overlap.b", "N]chr1:100]"),
        ],
    )

    assert [record["ID"] for record in records] == [
        "overlap.a",
        "overlap.b",
        "overlap.a",
    ]
    assert all(record["SVTYPE"] == "TRA" for record in records)
