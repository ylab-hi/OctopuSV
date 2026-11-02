"""correct: caller-aware SUPPORT, one record per BND junction, SvABA layout.

Records are taken from real VCaP caller output (Stage 16 audit) and trimmed.
"""
from __future__ import annotations

from pathlib import Path

import pytest
from typer.testing import CliRunner

pytest.importorskip("natsort")

from octopusv.cli.cli import app
from octopusv.utils.caller_evidence import detect_caller, resolve_support

runner = CliRunner()
CONTIGS = "".join(f"##contig=<ID={c},length=250000000>\n" for c in ("chr1", "chr2", "chr7", "chr8", "chr12", "chr22"))


def _vcf(source, records, samples=("S1",), extra_header=""):
    source_line = f"##source={source}\n" if source else ""
    return ("##fileformat=VCFv4.2\n" + source_line + CONTIGS + extra_header
            + "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(samples) + "\n"
            + "".join(line + "\n" for line in records))


def _correct(tmp_path, text, name="input.vcf"):
    source = tmp_path / name
    output = tmp_path / (name + ".svcf")
    source.write_text(text)
    result = runner.invoke(app, ["correct", "-i", str(source), "-o", str(output)])
    return result, output


def _records(output):
    rows = []
    for line in output.read_text().splitlines():
        if not line or line.startswith("#"):
            continue
        fields = line.split("\t")
        info = dict(item.split("=", 1) for item in fields[7].split(";") if "=" in item)
        rows.append(dict(CHROM=fields[0], POS=fields[1], ID=fields[2], ALT=fields[4],
                         samples=fields[9:], **info))
    return rows


def _support(tmp_path, source, records, **kwargs):
    result, output = _correct(tmp_path, _vcf(source, records, **kwargs))
    assert result.exit_code == 0, result.output
    return {row["ID"]: row["SUPPORT"] for row in _records(output)}


# ---------------------------------------------------------------- SUPPORT

def test_lumpy_uses_info_su_not_pe_plus_sr(tmp_path):
    rec = ("chr1\t934095\t1\tN\t<DEL>\t.\t.\tSVTYPE=DEL;SVLEN=-760;END=934855;IMPRECISE;"
           "SU=13;PE=13;SR=0\tGT:SU:PE:SR\t./.:13:13:0")
    assert _support(tmp_path, "LUMPY", [rec]) == {"1": "13"}


def test_delly_uses_sample_dv_plus_rv(tmp_path):
    rec = ("chr1\t789498\tDEL00000003\tA\t<DEL>\t225\tPASS\tIMPRECISE;SVTYPE=DEL;"
           "SVMETHOD=EMBL.DELLYv1.3.1;END=224012409;PE=8;MAPQ=25;CT=3to5\t"
           "GT:GL:GQ:FT:RCL:RC:RCR:RDCN:DR:DV:RR:RV\t0/1:-17,0,-43:10000:PASS:1:2:3:18:9:7:0:0")
    precise = ("chr1\t1062056\tDEL00000020\tCTCCCTCCCTTG\tC\t480\tPASS\tPRECISE;SVTYPE=DEL;"
               "SVMETHOD=EMBL.DELLYv1.3.1;END=1062067;PE=0;SR=8\t"
               "GT:GL:GQ:FT:RCL:RC:RCR:RDCN:DR:DV:RR:RV\t0/1:-1,0,-4:100:PASS:1:2:3:2:0:0:20:8")
    assert _support(tmp_path, None, [rec, precise]) == {"DEL00000003": "7", "DEL00000020": "8"}


def test_manta_adds_paired_and_split_alt_and_keeps_true_zero(tmp_path):
    records = [
        "chr1\t998763\tMantaINS:1\tGC\tGGGGAGGGCGCGGAGCGGAGGG\t.\tPASS\tEND=998764;SVTYPE=INS;SVLEN=20\tPR:SR\t0,0:0,32",
        "chr1\t1108845\tMantaDEL:2\tTCGGCACCATCCACCACGGCCACCATG\tT\t.\tPASS\tEND=1108871;SVTYPE=DEL;SVLEN=-26\tPR:SR\t5,0:14,6",
        "chr1\t2108845\tMantaDEL:3\tTCGGCACCATCCACCACGGCCACCATG\tT\t.\tPASS\tEND=2108871;SVTYPE=DEL;SVLEN=-26\tPR\t5,9",
        "chr1\t3108845\tMantaDEL:4\tTCGGCACCATCCACCACGGCCACCATG\tT\t.\tPASS\tEND=3108871;SVTYPE=DEL;SVLEN=-26\tPR\t5,0",
    ]
    assert _support(tmp_path, "GenerateSVCandidates 1.6.0", records) == {
        "MantaINS:1": "32", "MantaDEL:2": "6", "MantaDEL:3": "9", "MantaDEL:4": "0"}


def test_pbsv_uses_ad_alt(tmp_path):
    rec = ("chr1\t10832\tpbsv.DEL.1\tCGCAGGCGCAGAGAGGCGCGCCGCGCCGGC\tC\t.\tPASS\tSVTYPE=DEL;END=10861;SVLEN=-29\t"
           "GT:AD:DP:SAC\t0/1:23,14:37:14,9,7,7")
    assert _support(tmp_path, "pbsv 2.10.0 (commit v2.10.0)", [rec]) == {"pbsv.DEL.1": "14"}


def test_hitsv_long_read_alt_and_minus_one_is_unknown(tmp_path):
    records = [
        "chr1\t72128\th1\tT\tTATATATACACACATATATACATACATAC\t67\tPASS\tSVTYPE=INS;END=72129;SVLEN=28\t"
        "GT:GQ:LR:SR:SL:IF\t0/1:67:10,2,0:-1,-1,-1:28:11,0,2,0,1",
        "chr1\t82128\th2\tT\tTATATATACACACATATATACATACATAC\t67\tPASS\tSVTYPE=INS;END=82129;SVLEN=28\t"
        "GT:GQ:LR:SR:SL:IF\t0/1:67:-1,-1,-1:-1,-1,-1:28:11,0,2,0,1",
    ]
    assert _support(tmp_path, "HitSV2.0.4", records) == {"h1": "2", "h2": "."}


def test_unprofiled_callers_keep_historical_support(tmp_path):
    svim = ("chr1\t10862\tsvim.INS.3\tG\tGCAGGCGCAGAGAGGCGCGCCGCGCCGGCGCAGGCGCAG\t26\tPASS\t"
            "SVTYPE=INS;END=10862;SVLEN=38;SUPPORT=22\tGT:DP:AD\t0/1:28:6,22")
    assert _support(tmp_path, "SVIM-v2.0.0", [svim]) == {"svim.INS.3": "22"}
    pr_only = "chr1\t1108845\tx1\tTCGGCACCATCCACCACGGCCACCATG\tT\t.\tPASS\tEND=1108871;SVTYPE=DEL;SVLEN=-26\tPR:SR\t5,4:14,6"
    # Unknown caller: the historical first-sample PR rule is unchanged.
    assert _support(tmp_path, "SomeCaller", [pr_only]) == {"x1": "4"}


def test_multisample_support_sums_known_samples_and_keeps_sample_order(tmp_path):
    fmt = "GT:GL:GQ:FT:RCL:RC:RCR:RDCN:DR:DV:RR:RV"
    info = "PRECISE;SVTYPE=DEL;SVMETHOD=EMBL.DELLYv1.3.1;END=2000"
    both = f"chr1\t1000\td1\tN\t<DEL>\t60\tPASS\t{info}\t{fmt}\t0/1:0,0,0:9:PASS:1:2:3:2:9:4:9:3\t0/0:0,0,0:9:PASS:1:2:3:2:9:0:9:0"
    one_unknown = f"chr1\t5000\td2\tN\t<DEL>\t60\tPASS\t{info.replace('2000', '6000')}\t{fmt}\t1/1:0,0,0:9:PASS:1:2:3:2:9:4:9:3\t./.:.:.:.:.:.:.:.:.:.:.:."
    result, output = _correct(tmp_path, _vcf(None, [both, one_unknown], samples=("T", "N")))
    assert result.exit_code == 0, result.output
    rows = {row["ID"]: row for row in _records(output)}
    assert rows["d1"]["SUPPORT"] == "7"
    assert rows["d2"]["SUPPORT"] == "."
    header = [line for line in output.read_text().splitlines() if line.startswith("#CHROM")][0]
    assert header.split("\t")[9:] == ["T", "N"]
    assert rows["d1"]["samples"][0].startswith("0/1:") and rows["d1"]["samples"][1].startswith("0/0:")


def test_detect_and_resolve_units():
    assert detect_caller(".", {"SVMETHOD": "EMBL.DELLYv1.3.1"}) == "delly"
    assert detect_caller("GenerateSVCandidates", {}) == "manta"
    assert detect_caller("SVIM-v2.0.0", {}) is None
    assert resolve_support("lumpy", {"SU": "5"}, "GT:SU", ["./.:9"]) == "5"
    assert resolve_support("lumpy", {}, "GT:SU", ["./.:9", "./.:2"]) == "11"


# ---------------------------------------------------------------- BND junctions

def _bnd(chrom, pos, record_id, alt, info="SVTYPE=BND", fmt="GT", sample="0/1"):
    return f"{chrom}\t{pos}\t{record_id}\tN\t{alt}\t60\tPASS\t{info}\t{fmt}\t{sample}"


@pytest.mark.parametrize("alt_a,alt_b", [("N]chr2:200]", "N]chr1:100]"), ("[chr2:200[N", "[chr1:100[N"),
                                         ("N[chr2:200[", "]chr1:100]N")])
def test_two_ends_of_one_junction_give_one_record(tmp_path, alt_a, alt_b):
    result, output = _correct(tmp_path, _vcf("SVIM-v2.0.0", [
        _bnd("chr1", 100, "a", alt_a, "SVTYPE=BND;SUPPORT=5"),
        _bnd("chr2", 200, "b", alt_b, "SVTYPE=BND;SUPPORT=9"),
    ]))
    assert result.exit_code == 0, result.output
    rows = _records(output)
    assert len(rows) == 1 and rows[0]["SVTYPE"] == "TRA"
    assert rows[0]["SUPPORT"] == "9"  # same reads: larger value, never 5 + 9


@pytest.mark.parametrize("alt_a,alt_b", [("N[chr2:200[", "N[chr1:100["), ("]chr2:200]N", "]chr1:100]N"),
                                         ("N]chr2:200]", "[chr1:100[N"), ("N[chr2:200[", "N]chr1:100]"),
                                         ("N[chr2:200[", "[chr1:100[N"), ("]chr2:200]N", "N]chr1:100]"),
                                         ("]chr2:200]N", "[chr1:100[N")])
def test_two_different_junctions_keep_two_records(tmp_path, alt_a, alt_b):
    result, output = _correct(tmp_path, _vcf("SVIM-v2.0.0", [
        _bnd("chr1", 100, "a", alt_a), _bnd("chr2", 200, "b", alt_b)]))
    assert result.exit_code == 0, result.output
    rows = _records(output)
    assert sorted(row["ID"] for row in rows) == ["a", "b"]
    assert all(row["SVTYPE"] == "TRA" for row in rows)


def test_mateid_pairs_the_declared_mate_among_near_candidates(tmp_path):
    lumpy = lambda chrom, pos, rid, alt, mate: _bnd(
        chrom, pos, rid, alt, f"SVTYPE=BND;EVENT=1;MATEID={mate};SU=7;PE=7;SR=0", "GT:SU:PE:SR", "./.:7:7:0")
    result, output = _correct(tmp_path, _vcf("LUMPY", [
        lumpy("chr1", 100, "1_1", "N]chr2:200]", "1_2"),
        lumpy("chr2", 201, "2_1", "N]chr1:101]", "2_2"),   # near decoy, different event
        lumpy("chr2", 200, "1_2", "N]chr1:100]", "1_1"),
        lumpy("chr1", 101, "2_2", "N]chr2:201]", "2_1"),
    ]))
    assert result.exit_code == 0, result.output
    rows = _records(output)
    assert len(rows) == 2 and {row["SUPPORT"] for row in rows} == {"7"}
    assert "same junction:   2" in result.output


# ---------------------------------------------------------------- SvABA

SVABA_FMT = "GT:AD:DP:GQ:PL:SR:DR:LR:LO"


def _svaba(extra="0\t0\t12", sample="0/1:12:111:8.6:8.6,0,269.3:0:12:-8.411:23.5", pos=791524, rid="784459337:1",
           alt="]chr2:791725]C"):
    return (f"chr1\t{pos}\t{rid}\tC\t{alt}\t8\tPASS\tEVDNC=DSCRD;IMPRECISE;MAPQ=18;SVTYPE=BND\t"
            f"{SVABA_FMT}\t{extra}\t{sample}")


def test_svaba_known_layout_reads_last_column_and_keeps_extra_columns(tmp_path):
    text = _vcf("svaba(v1.1.0) run -t ../VCaP.bam", [_svaba()], samples=("../VCaP.bam",))
    result, output = _correct(tmp_path, text)
    assert result.exit_code == 0, result.output
    assert _records(output)[0]["SUPPORT"] == "12"
    side = Path(f"{output}.svaba_extra_columns.tsv").read_text().splitlines()
    assert side[1].split("\t") == ["10", "chr1", "791524", "784459337:1", "0", "0", "12"]
    assert "SvABA layout" in result.output


@pytest.mark.parametrize("extra,sample", [("0\tx\t12", "0/1:12:111:8.6:8.6,0,269.3:0:12:-8.411:23.5"),
                                          ("0\t0\t12", "0/1:12:111")])
def test_svaba_unverified_layout_is_rejected_and_writes_nothing(tmp_path, extra, sample):
    text = _vcf("svaba(v1.1.0)", [_svaba(extra=extra, sample=sample)], samples=("../VCaP.bam",))
    result, output = _correct(tmp_path, text)
    assert result.exit_code != 0
    assert not output.exists() and not Path(f"{output}.svaba_extra_columns.tsv").exists()


# ---------------------------------------------------------------- v2 additions

def _correct_args(tmp_path, text, args):
    source = tmp_path / "input.vcf"
    output = tmp_path / "input.vcf.svcf"
    source.write_text(text)
    return runner.invoke(app, ["correct", "-i", str(source), "-o", str(output), *args]), output


def test_unknown_caller_keeps_historical_support_filtering(tmp_path):
    rec = "chr1\t1000\tu1\tN\t<DEL>\t60\tPASS\tSVTYPE=DEL;END=2000;SVLEN=-1000\tGT:AD\t0/1:7,3"
    result, output = _correct_args(tmp_path, _vcf("SomeCaller", [rec]), ["--min-support", "10"])
    assert result.exit_code == 0, result.output
    assert _records(output) == []  # historical AD fallback (3) is below 10


def test_normal_multisample_svaba_keeps_every_sample(tmp_path):
    rec = ("chr1\t1000\ts1\tC\t]chr2:5000]C\t30\tPASS\tSVTYPE=BND\tGT:AD\t"
           "0/1:5\t0/0:0\t0/1:8\t0/1:13")
    result, output = _correct(tmp_path, _vcf("svaba(v1.1.0)", [rec], samples=("A", "B", "C", "D")))
    assert result.exit_code == 0, result.output
    row = _records(output)[0]
    assert row["SUPPORT"] == "26" and len(row["samples"]) == 4
    assert not Path(f"{output}.svaba_extra_columns.tsv").exists()


def test_mateid_with_contradicting_orientation_is_not_trusted(tmp_path):
    rec = lambda chrom, pos, rid, alt, mate: _bnd(
        chrom, pos, rid, alt, f"SVTYPE=BND;MATEID={mate};SU=4", "GT:SU", "./.:4")
    result, output = _correct(tmp_path, _vcf("LUMPY", [
        rec("chr1", 100, "m1", "N[chr2:200[", "m2"),
        rec("chr2", 200, "m2", "N[chr1:100[", "m1"),   # t[p[ + t[p[: two junctions
    ]))
    assert result.exit_code == 0, result.output
    assert sorted(row["ID"] for row in _records(output)) == ["m1", "m2"]
    assert "MATEID orientation conflicts:    1" in result.output


def test_multisample_mates_take_per_sample_max_then_sum(tmp_path):
    rec = lambda chrom, pos, rid, alt, s1, s2: (
        f"{chrom}\t{pos}\t{rid}\tN\t{alt}\t.\tPASS\tSVTYPE=BND\tGT:AD\t0/1:{s1}\t0/1:{s2}")
    result, output = _correct(tmp_path, _vcf("pbsv 2.10.0", [
        rec("chr1", 100, "p1", "N]chr2:200]", "5,10", "5,0"),
        rec("chr2", 200, "p2", "N]chr1:100]", "5,0", "5,10"),
    ], samples=("T", "N")))
    assert result.exit_code == 0, result.output
    rows = _records(output)
    assert len(rows) == 1 and rows[0]["SUPPORT"] == "20"


def test_issue_193_min_svlen_keeps_translocations_and_filters_linear_svs(tmp_path):
    records = [
        "chr1\t100\tdb.tra\tN\t<TRA>\t.\tPASS\tSVTYPE=TRA;CHR2=chr2;END=500;SVLEN=0;SUPPORT=8\tGT\t0/1",
        _bnd("chr1", 300, "bnd.single", "N[chr7:900[", "SVTYPE=BND;SUPPORT=6"),
        "chr1\t1000\tdel.small\tN\t<DEL>\t.\tPASS\tSVTYPE=DEL;END=1030;SVLEN=-30;SUPPORT=8\tGT\t0/1",
        "chr1\t2000\tdel.big\tN\t<DEL>\t.\tPASS\tSVTYPE=DEL;END=2100;SVLEN=-100;SUPPORT=8\tGT\t0/1",
        "chr1\t3000\tinv.big\tN\t<INV>\t.\tPASS\tSVTYPE=INV;END=3100;SVLEN=100;SUPPORT=8\tGT\t0/1",
        "chr1\t4000\tdup.big\tN\t<DUP>\t.\tPASS\tSVTYPE=DUP;END=4100;SVLEN=100;SUPPORT=8\tGT\t0/1",
    ]
    result, output = _correct_args(tmp_path, _vcf("DeBreak", records), ["--min-svlen", "50"])
    assert result.exit_code == 0, result.output
    rows = {row["ID"]: row for row in _records(output)}
    assert set(rows) == {"db.tra", "bnd.single", "del.big", "inv.big", "dup.big"}
    assert rows["db.tra"]["SVTYPE"] == "TRA"


def test_stale_svaba_side_file_is_removed(tmp_path):
    output = tmp_path / "input.vcf.svcf"
    stale = Path(f"{output}.svaba_extra_columns.tsv")
    stale.write_text("old\n")
    rec = "chr1\t1000\tx\tN\t<DEL>\t60\tPASS\tSVTYPE=DEL;END=2000;SVLEN=-1000;SUPPORT=4\tGT\t0/1"
    result, _ = _correct(tmp_path, _vcf("SVIM-v2.0.0", [rec]))
    assert result.exit_code == 0, result.output
    assert not stale.exists()
