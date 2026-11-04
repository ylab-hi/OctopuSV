"""Guard tests for sample-mode merge grouping (written before issue #192).

These lock the CURRENT sample-mode behavior for two same-sample situations:

1. One sample file holds both ends of one translocation junction as two
   records.  They describe one physical structure and must stay in one merged
   record, now and after #192.

2. One sample file holds two different deletions that each match another
   sample's seed but cannot be merged with each other directly (length ratio
   3.0).  Today they are merged into one record and one of them disappears
   from SOURCE_IDS.  Issue #192 is expected to split them; the test for that
   case documents the status quo and must be updated deliberately in the #192
   change, with the reason stated.
"""

from __future__ import annotations

from pathlib import Path

from typer.testing import CliRunner

from octopusv.cli.cli import app


runner = CliRunner()


def _invoke(args):
    result = runner.invoke(app, [str(arg) for arg in args])
    assert result.exit_code == 0, f"octopusv {' '.join(map(str, args))}\n{result.output}"
    return result


def _raw_vcf(path: Path, source: str, rows: list[str]) -> Path:
    path.write_text(
        "##fileformat=VCFv4.2\n"
        f"##source={source}\n"
        "##contig=<ID=chr1,length=1000000>\n"
        "##contig=<ID=chr2,length=1000000>\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n"
        + "".join(row + "\n" for row in rows)
    )
    return path


def _del(pos: int, record_id: str, length: int, gt: str = "0/1") -> str:
    return (
        f"chr1\t{pos}\t{record_id}\tN\t<DEL>\t60\tPASS\t"
        f"SVTYPE=DEL;END={pos + length};SVLEN=-{length};SUPPORT=8\tGT\t{gt}"
    )


def _bnd(chrom: str, pos: int, record_id: str, alt: str, gt: str = "0/1") -> str:
    return f"{chrom}\t{pos}\t{record_id}\tN\t{alt}\t60\tPASS\tSVTYPE=BND;SUPPORT=8\tGT\t{gt}"


def _sample_svcf(tmp_path: Path, name: str, rows: list[str]) -> Path:
    raw = _raw_vcf(tmp_path / f"{name}.vcf", "Sniffles2", rows)
    svcf = tmp_path / f"{name}.svcf"
    _invoke(["correct", "-i", raw, "-o", svcf])
    return svcf


def _merge_samples(tmp_path: Path, inputs: list[Path], names: str) -> list[list[str]]:
    out = tmp_path / "merged.svcf"
    args = ["merge"]
    for path in inputs:
        args += ["-i", path]
    args += ["--mode", "sample", "--sample-names", names, "--union", "-o", out]
    _invoke(args)
    _invoke(["validate-svcf", out])
    return [
        line.split("\t")
        for line in out.read_text().splitlines()
        if line and not line.startswith("#")
    ]


def _info(fields: list[str]) -> dict[str, str]:
    return dict(item.split("=", 1) for item in fields[7].split(";") if "=" in item)


def _gt(column: str) -> str:
    return column.split(":", 1)[0]


def test_same_sample_two_ends_of_one_translocation_stay_one_record(tmp_path):
    # Sample A reports both ends of one junction; the ends are 5-10 bp apart,
    # so `correct` keeps them as two TRA records in A's file.
    sample_a = _sample_svcf(
        tmp_path,
        "A",
        [
            _bnd("chr1", 300000, "A.end1", "N[chr2:5000["),
            _bnd("chr2", 5010, "A.end2", "]chr1:300005]N"),
        ],
    )
    sample_b = _sample_svcf(
        tmp_path, "B", [_bnd("chr1", 300002, "B.end1", "N[chr2:5003[")]
    )

    rows = _merge_samples(tmp_path, [sample_a, sample_b], "A,B")

    assert len(rows) == 1
    info = _info(rows[0])
    assert info["SVTYPE"] == "TRA"
    assert (rows[0][0], rows[0][1], rows[0][2]) == ("chr1", "300000", "A.end1")
    assert info["SOURCES"] == "A,B"
    assert info["SOURCE_IDS"] == "A.end1,B.end1"
    assert [_gt(column) for column in rows[0][9:]] == ["0/1", "0/1"]


def test_same_sample_two_different_deletions_status_quo_before_192(tmp_path):
    # B's 80 bp DEL is the seed (leftmost start).  A's 50 bp and 150 bp DELs
    # each match the seed (ratios 1.6 and 1.875) but not each other (3.0).
    sample_a = _sample_svcf(
        tmp_path,
        "A",
        [_del(10000, "A.del50", 50, "0/1"), _del(10000, "A.del150", 150, "1/1")],
    )
    sample_b = _sample_svcf(tmp_path, "B", [_del(9990, "B.del80", 80)])

    rows = _merge_samples(tmp_path, [sample_a, sample_b], "A,B")

    # STATUS QUO (issue #192): both A deletions fall into B's group; only one
    # A record ID survives and A's genotype is a consensus of two different
    # structures.  #192 must change this test on purpose: expect two records,
    # A.del50 and A.del150 kept apart and linked.
    assert len(rows) == 1
    info = _info(rows[0])
    assert info["SOURCES"] == "A,B"
    assert info["SOURCE_IDS"] == "A.del50,B.del80"
    assert "A.del150" not in rows[0][7]
    assert [_gt(column) for column in rows[0][9:]] == ["1/.", "0/1"]
