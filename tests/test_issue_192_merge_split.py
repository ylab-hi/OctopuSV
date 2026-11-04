"""Issue #192: records from one source that are different structures are kept
apart; records that are the same structure stay merged.

Caller mode rule for DEL/DUP/INS: two records from one caller are the same
structure when their lengths differ by at most max(10 bp, 10%).  INV, TRA and
BND groups are not split in caller mode.
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


def _caller_svcf(tmp_path: Path, name: str, source: str, rows: list[str]) -> Path:
    raw = tmp_path / f"{name}.vcf"
    raw.write_text(
        "##fileformat=VCFv4.2\n"
        f"##source={source}\n"
        "##contig=<ID=chr1,length=1000000>\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n"
        + "".join(row + "\n" for row in rows)
    )
    svcf = tmp_path / f"{name}.svcf"
    _invoke(["correct", "-i", raw, "-o", svcf])
    return svcf


def _sv(svtype: str, pos: int, record_id: str, length: int) -> str:
    end = pos if svtype == "INS" else pos + length
    svlen = -length if svtype == "DEL" else length
    return (
        f"chr1\t{pos}\t{record_id}\tN\t<{svtype}>\t60\tPASS\t"
        f"SVTYPE={svtype};END={end};SVLEN={svlen};SUPPORT=8\tGT\t0/1"
    )


def _merge(tmp_path: Path, inputs: list[Path], labels: str, name="merged") -> list[dict]:
    out = tmp_path / f"{name}.svcf"
    args = ["merge"]
    for path in inputs:
        args += ["-i", path]
    args += ["--mode", "caller", "--caller-names", labels, "--union", "-o", out]
    _invoke(args)
    _invoke(["validate-svcf", out])
    records = []
    for line in out.read_text().splitlines():
        if not line or line.startswith("#"):
            continue
        fields = line.split("\t")
        info = dict(item.split("=", 1) for item in fields[7].split(";") if "=" in item)
        records.append({"id": fields[2], "info": info,
                        "ids": frozenset(info["SOURCE_IDS"].split(","))})
    return records


def test_insertions_of_clearly_different_length_are_split_and_linked(tmp_path):
    a = _caller_svcf(tmp_path, "A", "cuteSV", [_sv("INS", 610547, "a.203", 203),
                                               _sv("INS", 610559, "a.71", 71)])
    b = _caller_svcf(tmp_path, "B", "SVIM", [_sv("INS", 610547, "b.204", 204),
                                             _sv("INS", 610591, "b.71", 71)])
    c = _caller_svcf(tmp_path, "C", "pbsv", [_sv("INS", 610550, "c.200", 200)])

    records = _merge(tmp_path, [a, b, c], "A,B,C")

    assert {record["ids"] for record in records} == {
        frozenset({"a.203", "b.204", "c.200"}),
        frozenset({"a.71", "b.71"}),
    }
    assert len({record["info"]["MERGE_SPLIT"] for record in records}) == 1


def test_same_caller_lengths_within_ten_percent_stay_merged(tmp_path):
    a = _caller_svcf(tmp_path, "A", "cuteSV", [_sv("DEL", 1000, "a.100", 100),
                                               _sv("DEL", 1004, "a.108", 108)])
    b = _caller_svcf(tmp_path, "B", "SVIM", [_sv("DEL", 1002, "b.104", 104)])

    records = _merge(tmp_path, [a, b], "A,B")

    assert len(records) == 1
    assert records[0]["ids"] == {"a.100", "a.108", "b.104"}
    assert "MERGE_SPLIT" not in records[0]["info"]


def test_small_length_difference_uses_the_10_bp_floor(tmp_path):
    # 50 vs 59 bp is 15% but only 9 bp: same structure.
    a = _caller_svcf(tmp_path, "A", "cuteSV", [_sv("DEL", 1000, "a.50", 50),
                                               _sv("DEL", 1000, "a.59", 59)])
    records = _merge(tmp_path, [a], "A")
    assert len(records) == 1 and "MERGE_SPLIT" not in records[0]["info"]


def test_inversion_halves_from_one_caller_are_not_split(tmp_path):
    a = _caller_svcf(tmp_path, "A", "DELLY", [_sv("INV", 5000, "a.half1", 812),
                                              _sv("INV", 5010, "a.half2", 700)])
    b = _caller_svcf(tmp_path, "B", "LUMPY", [_sv("INV", 5003, "b.inv", 805)])
    records = _merge(tmp_path, [a, b], "A,B")
    assert len(records) == 1
    assert records[0]["ids"] == {"a.half1", "a.half2", "b.inv"}


def test_split_membership_does_not_depend_on_input_order(tmp_path):
    a = _caller_svcf(tmp_path, "A", "cuteSV", [_sv("DEL", 1000, "a.56", 56),
                                               _sv("DEL", 1000, "a.112", 112)])
    b = _caller_svcf(tmp_path, "B", "SVIM", [_sv("DEL", 1001, "b.56", 56),
                                             _sv("DEL", 1001, "b.112", 112)])
    c = _caller_svcf(tmp_path, "C", "Manta", [_sv("DEL", 1000, "c.112", 112)])

    forward = _merge(tmp_path, [a, b, c], "A,B,C", "forward")
    reverse = _merge(tmp_path, [c, b, a], "C,B,A", "reverse")

    expected = {frozenset({"a.56", "b.56"}), frozenset({"a.112", "b.112", "c.112"})}
    assert {record["ids"] for record in forward} == expected
    assert {record["ids"] for record in reverse} == expected


def test_record_between_two_split_structures_is_not_forced_into_either(tmp_path):
    # A reports 100 and 200 bp; B's 150 bp is equivalent to neither, so it is
    # kept as its own record rather than attached to the nearer one.  C's
    # 152 bp is equivalent to B's and joins it.
    a = _caller_svcf(tmp_path, "A", "cuteSV", [_sv("DEL", 1000, "a.100", 100),
                                               _sv("DEL", 1000, "a.200", 200)])
    b = _caller_svcf(tmp_path, "B", "SVIM", [_sv("DEL", 1000, "b.150", 150)])
    c = _caller_svcf(tmp_path, "C", "pbsv", [_sv("DEL", 1000, "c.152", 152)])

    records = _merge(tmp_path, [a, b, c], "A,B,C")

    assert {record["ids"] for record in records} == {
        frozenset({"a.100"}),
        frozenset({"a.200"}),
        frozenset({"b.150", "c.152"}),
    }
    assert len({record["info"]["MERGE_SPLIT"] for record in records}) == 1
