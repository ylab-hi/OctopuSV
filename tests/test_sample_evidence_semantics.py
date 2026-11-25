"""Sample-mode "has evidence" is one fact, read the same way by every command.

In SVCF 1.1 sample mode a column has evidence unless it is the writer's
no-event placeholder (GT=0/0, UC=0, UV=0, no ID/SC).  A caller's reference
call (GT=0/0, e.g. SVDSS) and an unresolved call (GT=./., e.g. SVIM) are
evidence.  subset used to drop such records and strip them from SOURCES;
inspect flagged them as placeholders and warned about a false mismatch.
"""

import json

import pytest
from typer.testing import CliRunner

from octopusv.cli.cli import app
from octopusv.utils.sample_mode_semantics import sample_column_has_evidence

FORMAT11 = "GT:AD:UC:UV:LN:ST:QV:TY:ID:SC:REF:ALT:CO"
LEGACY = "GT:AD:LN:ST:QV:TY:ID:SC:REF:ALT:CO"
TAIL = "500:.:30:DEL:{id}:{sc}:N:<DEL>:chr1_1000-chr1_1500"


def _block(gt, uc, uv, record_id="A.1", sc="OctopuSV"):
    return f"{gt}:.,.:{uc}:{uv}:" + TAIL.format(id=record_id, sc=sc)


@pytest.mark.parametrize(
    "block, expected",
    [
        (_block("0/1", 1, 1), True),                      # carrier
        (_block("0/0", 0, 1), True),                      # caller's reference call
        (_block("./.", 0, 0), True),                      # evidence, no valid vote
        (_block("0/0", 0, 0, record_id=".", sc="."), False),  # writer placeholder
        ("0/0:.:.:.:.:.:.:.:.:.:.:.:.", False),           # column padding
        (None, False),
        (".", False),
    ],
)
def test_svcf11_sample_column_evidence(block, expected):
    assert sample_column_has_evidence(FORMAT11, block) is expected


def test_legacy_layout_keeps_gt_reading():
    block = "0/0:.:500:.:30:DEL:A.1:sniffles:N:<DEL>:chr1_1000-chr1_1500"
    assert sample_column_has_evidence(LEGACY, block) is False
    assert sample_column_has_evidence(LEGACY, "0/1" + block[3:]) is True


def _single(path, sample, gt):
    path.write_text(
        "##fileformat=SVCFv1.1\n##contig=<ID=chr1,length=1000000>\n"
        f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample}\n"
        f"chr1\t1000\t{sample}.1\tN\t<DEL>\t30\tPASS\t"
        "SVTYPE=DEL;END=1500;SVLEN=500;CHR2=chr1;SUPPORT=10;SVMETHOD=OctopuSV;"
        "RTID=.;AF=.;STRAND=+-;RNAMES=.\t"
        f"{LEGACY}\t{gt}:5,10:500:+-:30:DEL:{sample}.1:caller{sample}:N:<DEL>:chr1_1000-chr1_1500\n"
    )


@pytest.fixture
def merged(tmp_path):
    """One merged DEL: sA has GT 0/0, sB GT ./., sC GT 0/1 -- all with evidence."""
    inputs = []
    for sample, gt in (("sA", "0/0"), ("sB", "./."), ("sC", "0/1")):
        path = tmp_path / f"{sample}.svcf"
        _single(path, sample, gt)
        inputs += ["-i", str(path)]
    out = tmp_path / "merged.svcf"
    result = CliRunner().invoke(app, ["merge", *inputs, "-o", str(out), "--mode", "sample", "--union"])
    assert result.exit_code == 0, result.output
    return out


def _records(path):
    return [line.rstrip("\n").split("\t") for line in path.read_text().splitlines()
            if not line.startswith("#")]


def _info(cols):
    return dict(item.split("=", 1) for item in cols[7].split(";") if "=" in item)


@pytest.mark.parametrize("sample, source_id", [("sA", "sA.1"), ("sB", "sB.1"), ("sC", "sC.1")])
def test_subset_keeps_evidence_backed_non_carrier(merged, tmp_path, sample, source_id):
    out = tmp_path / f"subset_{sample}.svcf"
    result = CliRunner().invoke(app, ["subset", "-i", str(merged), "-o", str(out), "--sample", sample])
    assert result.exit_code == 0, result.output

    rows = _records(out)
    assert len(rows) == 1
    info = _info(rows[0])
    assert info["SOURCES"] == sample
    assert info["SOURCE_IDS"] == source_id


def test_subset_all_samples_keeps_provenance(merged, tmp_path):
    out = tmp_path / "subset_all.svcf"
    args = ["subset", "-i", str(merged), "-o", str(out),
            "--sample", "sA", "--sample", "sB", "--sample", "sC"]
    result = CliRunner().invoke(app, args)
    assert result.exit_code == 0, result.output
    info = _info(_records(out)[0])
    assert info["SOURCES"] == "sA,sB,sC"
    assert info["SOURCE_IDS"] == "sA.1,sB.1,sC.1"


def test_subset_still_drops_unobserved_sample(tmp_path):
    # sA carries a DEL at 1000; sB has only an unrelated DEL far away, so in
    # the first record sB's column is the writer placeholder.
    _single(tmp_path / "sA.svcf", "sA", "0/1")
    text = (tmp_path / "sA.svcf").read_text().replace("sA", "sB")
    text = text.replace("chr1\t1000\t", "chr1\t500000\t").replace("END=1500", "END=500500")
    text = text.replace("chr1_1000-chr1_1500", "chr1_500000-chr1_500500")
    (tmp_path / "sB.svcf").write_text(text)
    merged = tmp_path / "m.svcf"
    result = CliRunner().invoke(app, ["merge", "-i", str(tmp_path / "sA.svcf"), "-i", str(tmp_path / "sB.svcf"),
                                      "-o", str(merged), "--mode", "sample", "--union"])
    assert result.exit_code == 0, result.output
    assert len(_records(merged)) == 2

    out = tmp_path / "only_sB.svcf"
    result = CliRunner().invoke(app, ["subset", "-i", str(merged), "-o", str(out), "--sample", "sB"])
    assert result.exit_code == 0, result.output
    rows = _records(out)
    assert len(rows) == 1
    assert rows[0][1] == "500000"


def test_inspect_flags_evidence_not_genotype(merged):
    result = CliRunner().invoke(app, ["inspect", "-i", str(merged), "--id", "sA.1", "--json"])
    assert result.exit_code == 0, result.output
    record = json.loads(result.output)["results"][0]["record"]

    flags = {e["sample"]: (e["has_evidence"], e["is_placeholder"], e["is_carrier"])
             for e in record["evidence"]}
    assert flags == {
        "sA": (True, False, False),
        "sB": (True, False, False),
        "sC": (True, False, True),
    }
    assert record["evidence_samples"] == ["sA", "sB", "sC"]
    assert record["carrier_samples"] == ["sC"]
    assert not any(w["code"].startswith("W_SAMPLE_SOURCES_") for w in record["warnings"])


@pytest.mark.parametrize(
    "gt, carrier",
    [("0/0", False), ("0|0", False), ("0", False), ("./.", False), (".|.", False),
     (".", False), ("0/2", False), ("0/1", True), ("1|0", True), ("1/1", True),
     ("1", True), ("1/.", True)],
)
def test_genotype_carries_alt_uses_consensus_rules(gt, carrier):
    from octopusv.utils.sample_consensus import genotype_carries_alt

    assert genotype_carries_alt(gt) is carrier


@pytest.mark.parametrize("gt", ["0|0", "0", ".|."])
def test_inspect_evidence_backed_non_carrier_gts(merged, tmp_path, gt):
    # Rewrite sA's evidence-backed 0/0 as another non-carrier GT spelling.
    lines = merged.read_text().splitlines()
    for i, line in enumerate(lines):
        if not line.startswith("#"):
            cols = line.split("\t")
            assert cols[9].startswith("0/0:")
            cols[9] = gt + cols[9][3:]
            lines[i] = "\t".join(cols)
    edited = tmp_path / f"edited.svcf"
    edited.write_text("\n".join(lines) + "\n")

    result = CliRunner().invoke(app, ["inspect", "-i", str(edited), "--id", "sA.1", "--json"])
    assert result.exit_code == 0, result.output
    record = json.loads(result.output)["results"][0]["record"]
    sample_a = record["evidence"][0]
    assert sample_a["has_evidence"] is True
    assert sample_a["is_carrier"] is False
    assert record["carrier_samples"] == ["sC"]
