from __future__ import annotations

from pathlib import Path

import pytest
from typer.testing import CliRunner

pytest.importorskip("natsort")

from octopusv.cli.cli import app
from octopusv.utils.svcf_sample_parser import parse_svcf_sample_block


runner = CliRunner()


def _invoke(args: list[str]):
    result = runner.invoke(app, args)
    assert result.exit_code == 0, (
        f"Command failed: octopusv {' '.join(args)}\n"
        f"output:\n{result.output}\n"
        f"exception: {result.exception!r}"
    )
    return result


def _write_raw_vcf(
    path: Path,
    *,
    source: str,
    records: list[str],
    sample: str = "S1",
):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        "".join(
            [
                "##fileformat=VCFv4.2\n",
                f"##source={source}\n",
                "##contig=<ID=chr1,length=1000000>\n",
                "##contig=<ID=chr2,length=1000000>\n",
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"
                f"{sample}\n",
                *records,
            ]
        ),
        encoding="utf-8",
    )


def _raw_sv(
    chrom: str,
    pos: int,
    record_id: str,
    svtype: str,
    *,
    end: int,
    svlen: int | str,
    gt: str = "0/1",
    chr2: str | None = None,
    alt: str | None = None,
    support: int = 10,
):
    if alt is None:
        alt = f"<{svtype}>"
    info = [
        f"SVTYPE={svtype}",
        f"END={end}",
        f"SVLEN={svlen}",
        f"RE={support}",
    ]
    if chr2 is not None:
        info.append(f"CHR2={chr2}")
    return (
        f"{chrom}\t{pos}\t{record_id}\tN\t{alt}\t60\tPASS\t"
        f"{';'.join(info)}\tGT\t{gt}\n"
    )


def _raw_bnd(
    chrom: str,
    pos: int,
    record_id: str,
    alt: str,
    *,
    gt: str = "0/1",
):
    return (
        f"{chrom}\t{pos}\t{record_id}\tN\t{alt}\t60\tPASS\t"
        f"SVTYPE=BND;RE=10\tGT\t{gt}\n"
    )


def _correct(raw_vcf: Path, output_svcf: Path):
    output_svcf.parent.mkdir(parents=True, exist_ok=True)
    _invoke(["correct", "-i", str(raw_vcf), "-o", str(output_svcf)])
    _invoke(["validate-svcf", str(output_svcf)])
    return output_svcf


def _records(path: Path):
    return [
        line.split("\t")
        for line in path.read_text(encoding="utf-8").splitlines()
        if line and not line.startswith("#")
    ]


def _info(text: str):
    parsed = {}
    for item in text.split(";"):
        if "=" in item:
            key, value = item.split("=", 1)
            parsed[key] = value
        else:
            parsed[item] = True
    return parsed


def _by_id(path: Path):
    return {row[2]: row for row in _records(path)}


def _vcf_sample_map(row: list[str]):
    keys = row[8].split(":")
    return [dict(zip(keys, sample.split(":"), strict=False)) for sample in row[9:]]


def _find_record_near(path: Path, *, svtype: str, pos_lo: int, pos_hi: int):
    matches = []
    for row in _records(path):
        info = _info(row[7])
        if info.get("SVTYPE") == svtype and pos_lo <= int(row[1]) <= pos_hi:
            matches.append(row)
    assert len(matches) == 1, [row[:8] for row in matches]
    return matches[0]


def test_real_raw_to_final_sample_pipeline_preserves_binding_and_missingness(tmp_path):
    """Lock the real user path from raw caller VCF through final cohort VCF."""
    raw_a1 = tmp_path / "raw" / "sampleA.caller1.vcf"
    raw_a2 = tmp_path / "raw" / "sampleA.caller2.vcf"
    raw_b1 = tmp_path / "raw" / "sampleB.caller1.vcf"
    raw_b2 = tmp_path / "raw" / "sampleB.caller2.vcf"

    _write_raw_vcf(
        raw_a1,
        source="callerA1",
        records=[
            _raw_sv("chr1", 100, "A1.common", "DEL", end=600, svlen=-500, gt="0/1"),
            _raw_sv("chr1", 5000, "A1.only", "DEL", end=5500, svlen=-500, gt="0/1"),
        ],
    )
    _write_raw_vcf(
        raw_a2,
        source="callerA2",
        records=[
            _raw_sv("chr1", 105, "A2.common", "DEL", end=605, svlen=-500, gt="1/1"),
        ],
    )
    _write_raw_vcf(
        raw_b1,
        source="callerB1",
        records=[
            _raw_sv("chr1", 101, "B1.common", "DEL", end=601, svlen=-500, gt="0/1"),
        ],
    )
    _write_raw_vcf(
        raw_b2,
        source="callerB2",
        records=[
            _raw_sv("chr1", 103, "B2.common", "DEL", end=603, svlen=-500, gt="0/1"),
        ],
    )

    a1 = _correct(raw_a1, tmp_path / "corrected" / "A1.svcf")
    a2 = _correct(raw_a2, tmp_path / "corrected" / "A2.svcf")
    b1 = _correct(raw_b1, tmp_path / "corrected" / "B1.svcf")
    b2 = _correct(raw_b2, tmp_path / "corrected" / "B2.svcf")

    sample_a = tmp_path / "sampleA.caller.svcf"
    _invoke(
        [
            "merge",
            "-i", str(a1),
            "-i", str(a2),
            "-o", str(sample_a),
            "--mode", "caller",
            "--caller-names", "callerA1,callerA2",
            "--union",
        ]
    )

    sample_b = tmp_path / "sampleB.caller.svcf"
    _invoke(
        [
            "merge",
            "-i", str(b1),
            "-i", str(b2),
            "-o", str(sample_b),
            "--mode", "caller",
            "--caller-names", "callerB1,callerB2",
            "--union",
        ]
    )

    # Real caller-merge output must preserve exact source/evidence binding.
    common_a = _find_record_near(sample_a, svtype="DEL", pos_lo=95, pos_hi=110)
    common_a_info = _info(common_a[7])
    assert common_a_info["SOURCES"] == "callerA1,callerA2"
    assert common_a_info["SOURCE_IDS"] == "A1.common,A2.common"
    blocks = [
        parse_svcf_sample_block(common_a[8], block)
        for block in common_a[9:]
    ]
    assert [block["ID"] for block in blocks] == ["A1.common", "A2.common"]
    assert [block["GT"] for block in blocks] == ["0/1", "1/1"]

    population = tmp_path / "population.svcf"
    _invoke(
        [
            "merge",
            "-i", str(sample_a),
            "-i", str(sample_b),
            "-o", str(population),
            "--mode", "sample",
            "--sample-names", "sampleA,sampleB",
            "--union",
        ]
    )

    validation = _invoke(["validate-svcf", str(population)])
    assert "Status: PASS" in validation.output

    common = _find_record_near(population, svtype="DEL", pos_lo=95, pos_hi=110)
    assert common[8] == "GT:AD:UC:UV:LN:ST:QV:TY:ID:SC:REF:ALT:CO"
    sample_a_call = parse_svcf_sample_block(common[8], common[9])
    sample_b_call = parse_svcf_sample_block(common[8], common[10])
    assert (sample_a_call["GT"], sample_a_call["UC"], sample_a_call["UV"]) == (
        "1/.", "2", "2"
    )
    assert (sample_b_call["GT"], sample_b_call["UC"], sample_b_call["UV"]) == (
        "0/1", "2", "2"
    )

    only_a = _find_record_near(population, svtype="DEL", pos_lo=4995, pos_hi=5005)
    only_a_a = parse_svcf_sample_block(only_a[8], only_a[9])
    only_a_b = parse_svcf_sample_block(only_a[8], only_a[10])
    assert (only_a_a["GT"], only_a_a["UC"], only_a_a["UV"]) == (
        "0/1", "1", "1"
    )
    assert (only_a_b["GT"], only_a_b["UC"], only_a_b["UV"]) == (
        "0/0", "0", "0"
    )
    assert only_a_b["ID"] == "."
    assert only_a_b["SC"] == "."

    final_vcf = tmp_path / "population.vcf"
    _invoke(["svcf2vcf", "-i", str(population), "-o", str(final_vcf)])

    final_common = _find_record_near(final_vcf, svtype="DEL", pos_lo=95, pos_hi=110)
    final_common_samples = _vcf_sample_map(final_common)
    assert (final_common_samples[0]["GT"], final_common_samples[0]["UC"], final_common_samples[0]["UV"]) == (
        "1/.", "2", "2"
    )
    assert (final_common_samples[1]["GT"], final_common_samples[1]["UC"], final_common_samples[1]["UV"]) == (
        "0/1", "2", "2"
    )

    final_only_a = _find_record_near(final_vcf, svtype="DEL", pos_lo=4995, pos_hi=5005)
    final_only_a_samples = _vcf_sample_map(final_only_a)
    assert final_only_a_samples[0]["GT"] == "0/1"
    # Default VCF export must not turn an unobserved sample placeholder into an
    # evidence-backed homozygous-reference claim.
    assert final_only_a_samples[1]["GT"] == "./."
    assert final_only_a_samples[1]["UC"] == "0"
    assert final_only_a_samples[1]["UV"] == "0"


def test_same_source_multiple_evidence_counts_once_for_support_and_consensus(tmp_path):
    """Evidence multiplicity must not become caller support or extra genotype votes."""
    raw_a = tmp_path / "rawA.vcf"
    raw_b = tmp_path / "rawB.vcf"
    raw_c = tmp_path / "rawC.vcf"

    _write_raw_vcf(
        raw_a,
        source="sourceA",
        records=[
            _raw_sv("chr1", 100, "A.het", "DEL", end=600, svlen=-500, gt="0/1"),
            _raw_sv("chr1", 105, "A.hom", "DEL", end=605, svlen=-500, gt="1/1"),
        ],
    )
    _write_raw_vcf(
        raw_b,
        source="sourceB",
        records=[
            _raw_sv("chr1", 102, "B.het", "DEL", end=602, svlen=-500, gt="0/1"),
        ],
    )
    _write_raw_vcf(
        raw_c,
        source="sourceC",
        records=[
            _raw_sv("chr1", 20_000, "C.other", "DEL", end=20_500, svlen=-500, gt="0/1"),
        ],
    )

    a = _correct(raw_a, tmp_path / "A.svcf")
    b = _correct(raw_b, tmp_path / "B.svcf")
    c = _correct(raw_c, tmp_path / "C.svcf")

    min2 = tmp_path / "min2.svcf"
    _invoke(
        [
            "merge",
            "-i", str(a),
            "-i", str(b),
            "-i", str(c),
            "-o", str(min2),
            "--mode", "caller",
            "--caller-names", "sourceA,sourceB,sourceC",
            "--min-support", "2",
        ]
    )

    rows = _records(min2)
    assert len(rows) == 1
    merged = rows[0]
    merged_info = _info(merged[7])
    assert merged_info["SOURCES"] == "sourceA,sourceA,sourceB"
    assert merged_info["SOURCE_IDS"] == "A.het,A.hom,B.het"
    assert len(merged[9:]) == 3  # all evidence is preserved

    min3 = tmp_path / "min3.svcf"
    _invoke(
        [
            "merge",
            "-i", str(a),
            "-i", str(b),
            "-i", str(c),
            "-o", str(min3),
            "--mode", "caller",
            "--caller-names", "sourceA,sourceB,sourceC",
            "--min-support", "3",
        ]
    )
    # Three evidence blocks are only two independent supporting sources.
    assert _records(min3) == []

    final_vcf = tmp_path / "min2.vcf"
    _invoke(["svcf2vcf", "-i", str(min2), "-o", str(final_vcf)])
    final_row = _records(final_vcf)[0]
    call = _vcf_sample_map(final_row)[0]
    assert call["GT"] == "1/."
    assert call["UC"] == "2"
    assert call["UV"] == "2"
    assert call["AD"] == ".,."
    assert call["DP"] == "."


def test_all_svtypes_keep_geometry_semantics_across_correct_merge_and_vcf_export(tmp_path):
    """Lock representation changes without requiring byte-identical files."""
    raw = tmp_path / "all_types.vcf"
    _write_raw_vcf(
        raw,
        source="matrixCaller",
        records=[
            _raw_sv("chr1", 100, "native.del", "DEL", end=200, svlen=-100),
            _raw_sv("chr1", 1000, "native.dup", "DUP", end=1150, svlen=150),
            _raw_sv("chr1", 2000, "native.inv", "INV", end=2200, svlen=200),
            _raw_sv("chr1", 3000, "native.ins", "INS", end=3000, svlen=47),
            # Intentional design: one explicit interchromosomal BND becomes bracket TRA.
            _raw_bnd("chr1", 4000, "single.tra", "N[chr2:5000["),
            # Native symbolic TRA carries its remote breakpoint in CHR2 + END.
            _raw_sv(
                "chr1", 4500, "symbolic.tra", "TRA",
                end=5500, svlen=".", chr2="chr2", alt="<TRA>",
            ),
            # Same-chromosome explicit mate that cannot be safely reclassified stays BND.
            _raw_bnd("chr1", 6000, "kept.bnd", "N[chr1:6300["),
        ],
    )

    corrected = _correct(raw, tmp_path / "all_types.svcf")
    corrected_rows = _by_id(corrected)

    assert _info(corrected_rows["native.del"][7])["SVTYPE"] == "DEL"
    assert _info(corrected_rows["native.dup"][7])["SVTYPE"] == "DUP"
    assert _info(corrected_rows["native.inv"][7])["SVTYPE"] == "INV"

    ins_info = _info(corrected_rows["native.ins"][7])
    assert ins_info["SVTYPE"] == "INS"
    assert corrected_rows["native.ins"][1] == "3000"
    assert ins_info["SVLEN"] == "47"
    assert ins_info["END"] == "3047"  # SVCF 1.1 internal INS convention

    tra_info = _info(corrected_rows["single.tra"][7])
    assert tra_info["SVTYPE"] == "TRA"
    assert tra_info["CHR2"] == "chr2"
    assert tra_info["END"] == "5000"
    assert corrected_rows["single.tra"][4] == "N[chr2:5000["

    symbolic_tra_info = _info(corrected_rows["symbolic.tra"][7])
    assert symbolic_tra_info["SVTYPE"] == "TRA"
    assert symbolic_tra_info["CHR2"] == "chr2"
    assert symbolic_tra_info["END"] == "5500"
    assert corrected_rows["symbolic.tra"][4] == "<TRA>"

    bnd_info = _info(corrected_rows["kept.bnd"][7])
    assert bnd_info["SVTYPE"] == "BND"
    assert bnd_info["CHR2"] == "chr1"
    assert bnd_info["END"] == "6300"
    assert corrected_rows["kept.bnd"][4] == "N[chr1:6300["

    merged = tmp_path / "all_types.merged.svcf"
    _invoke(
        [
            "merge",
            "-i", str(corrected),
            "-o", str(merged),
            "--mode", "caller",
            "--caller-names", "matrixCaller",
            "--union",
        ]
    )
    validation = _invoke(["validate-svcf", str(merged)])
    assert "Status: PASS" in validation.output

    merged_rows = _by_id(merged)
    assert set(merged_rows) == set(corrected_rows)

    # Merge may choose/serialize representatives, but with one source record per
    # event the scientific geometry must be unchanged.
    for record_id in [
        "native.del", "native.dup", "native.inv",
        "single.tra", "symbolic.tra", "kept.bnd",
    ]:
        before = corrected_rows[record_id]
        after = merged_rows[record_id]
        before_info = _info(before[7])
        after_info = _info(after[7])
        assert (after[0], after[1], after[4]) == (before[0], before[1], before[4])
        assert after_info["SVTYPE"] == before_info["SVTYPE"]
        assert after_info["END"] == before_info["END"]
        assert after_info["CHR2"] == before_info["CHR2"]
        assert after_info["SVLEN"] == before_info["SVLEN"]

    merged_ins = _info(merged_rows["native.ins"][7])
    assert merged_ins["END"] == "3047"
    assert merged_ins["SVLEN"] == "47"

    final_vcf = tmp_path / "all_types.vcf.out"
    _invoke(["svcf2vcf", "-i", str(merged), "-o", str(final_vcf)])
    final_rows = _by_id(final_vcf)
    assert set(final_rows) == set(merged_rows)

    for record_id in ["native.del", "native.dup", "native.inv"]:
        before = merged_rows[record_id]
        after = final_rows[record_id]
        before_info = _info(before[7])
        after_info = _info(after[7])
        assert (after[0], after[1]) == (before[0], before[1])
        assert after_info["SVTYPE"] == before_info["SVTYPE"]
        assert after_info["END"] == before_info["END"]
        # SVCF stores positive absolute length. VCF export restores the
        # conventional negative sign for deletions, so compare event magnitude
        # rather than requiring byte-identical SVLEN text across formats.
        assert abs(int(after_info["SVLEN"])) == abs(int(before_info["SVLEN"]))
        if record_id == "native.del":
            assert int(after_info["SVLEN"]) < 0
        else:
            assert int(after_info["SVLEN"]) > 0

    final_ins_info = _info(final_rows["native.ins"][7])
    assert final_rows["native.ins"][1] == "3000"
    assert final_ins_info["SVTYPE"] == "INS"
    assert final_ins_info["SVLEN"] == "47"
    assert final_ins_info["END"] == "3000"  # conventional VCF INS endpoint

    final_tra = final_rows["single.tra"]
    final_tra_info = _info(final_tra[7])
    assert final_tra_info["SVTYPE"] == "TRA"
    assert final_tra_info["CHR2"] == "chr2"
    # Breakend-form TRA keeps the remote coordinate in ALT; INFO/END is not
    # required in conventional VCF export when ALT already encodes the mate.
    assert final_tra[4] == "N[chr2:5000["

    final_symbolic_tra = final_rows["symbolic.tra"]
    final_symbolic_tra_info = _info(final_symbolic_tra[7])
    assert final_symbolic_tra_info["SVTYPE"] == "TRA"
    assert final_symbolic_tra_info["CHR2"] == "chr2"
    assert final_symbolic_tra_info["END"] == "5500"
    assert final_symbolic_tra[4] == "<TRA>"

    final_bnd = final_rows["kept.bnd"]
    final_bnd_info = _info(final_bnd[7])
    assert final_bnd_info["SVTYPE"] == "BND"
    assert final_bnd_info["CHR2"] == "chr1"
    # BND mate position remains encoded in breakend ALT.
    assert final_bnd[4] == "N[chr1:6300["
