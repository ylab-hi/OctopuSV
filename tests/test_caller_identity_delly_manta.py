"""Caller identity that the input already states must reach SC.

DELLY writes no ##source line, but every record states SVMETHOD=EMBL.DELLY...
`correct` must write that value as FORMAT/SC (instead of ".") so the record
has an explicit caller identity: sample merge and source operations need it.
An explicit ##source
is never overridden and an unknown caller is never guessed.

Representative scoring is unchanged: DELLY was never ranked (its SC was
"."), and it stays unranked until the scoring itself is reviewed.
"""

from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace

from typer.testing import CliRunner

from octopusv.cli.cli import app
from octopusv.merger.sv_selector import select_representative_sv
from octopusv.utils.svcf_sample_parser import parse_svcf_sample_block


RUNNER = CliRunner()
DELLY_METHOD = "EMBL.DELLYv1.3.1"


def _delly_vcf(path: Path, source_line: str | None = "", method: str = DELLY_METHOD) -> Path:
    header = "##fileformat=VCFv4.2\n"
    if source_line:
        header += source_line + "\n"
    header += (
        "##contig=<ID=chr1,length=1000000>\n"
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tVCaP\n"
    )
    rows = [
        f"chr1\t10000\tDEL00000001\tN\t<DEL>\t1000\tPASS\tIMPRECISE;SVTYPE=DEL;SVMETHOD={method};"
        "END=12000;PE=10;MAPQ=60;CT=3to5;CIPOS=-50,50;CIEND=-50,50\tGT:GQ:DR:DV:RR:RV\t0/1:99:20:10:0:0",
        f"chr1\t50000\tDUP00000002\tN\t<DUP>\t800\tPASS\tIMPRECISE;SVTYPE=DUP;SVMETHOD={method};"
        "END=58000;PE=8;MAPQ=60;CT=5to3;CIPOS=-30,30;CIEND=-30,30\tGT:GQ:DR:DV:RR:RV\t0/1:90:18:8:0:0",
    ]
    path.write_text(header + "".join(row + "\n" for row in rows))
    return path


def _correct(tmp_path: Path, raw: Path) -> Path:
    out = tmp_path / (raw.stem + ".svcf")
    result = RUNNER.invoke(app, ["correct", "-i", str(raw), "-o", str(out)])
    assert result.exit_code == 0, result.output
    return out


def _sc_values(svcf: Path) -> list[str]:
    values = []
    for line in svcf.read_text().splitlines():
        if line and not line.startswith("#"):
            fields = line.split("\t")
            values.append(parse_svcf_sample_block(fields[8], fields[9])["SC"])
    return values


def test_delly_without_source_line_gets_its_svmethod_as_sc(tmp_path):
    out = _correct(tmp_path, _delly_vcf(tmp_path / "delly.vcf", source_line=None))
    assert _sc_values(out) == [DELLY_METHOD, DELLY_METHOD]


def test_delly_correct_output_now_works_in_sample_merge(tmp_path):
    (tmp_path / "a").mkdir()
    (tmp_path / "b").mkdir()
    a = _correct(tmp_path / "a", _delly_vcf(tmp_path / "a" / "s1.vcf", source_line=None))
    b = _correct(tmp_path / "b", _delly_vcf(tmp_path / "b" / "s2.vcf", source_line=None))
    merged = tmp_path / "merged.svcf"
    result = RUNNER.invoke(app, ["merge", "--mode", "sample", "--union", "-i", str(a), "-i", str(b),
                                 "-o", str(merged)])
    assert result.exit_code == 0, result.output
    assert RUNNER.invoke(app, ["validate-svcf", "-i", str(merged)]).exit_code == 0


def test_explicit_source_line_is_never_overridden(tmp_path):
    out = _correct(tmp_path, _delly_vcf(tmp_path / "delly.vcf", source_line="##source=MyDellyBuild"))
    assert _sc_values(out) == ["MyDellyBuild", "MyDellyBuild"]


def test_unknown_caller_without_source_line_is_not_guessed(tmp_path):
    out = _correct(tmp_path, _delly_vcf(tmp_path / "other.vcf", source_line=None, method="SomeCaller1.0"))
    assert _sc_values(out) == [".", "."]


def _member(caller: str, record_id: str, support: str = "10", quality: str = "50"):
    return SimpleNamespace(
        info={"SUPPORT": support}, quality=quality, sample={"SC": caller}, sv_id=record_id,
        source_file=f"/data/{record_id}.svcf", sample_name="S", format="GT",
    )


def test_named_delly_scores_exactly_as_unnamed_delly():
    # Representative scoring is unchanged: DELLY was unranked (SC=".") and
    # stays unranked now that it carries its name.
    for group in (
        # CHD1 truth (VCaP): DeBreak support 6, no QUAL; DELLY support 2, QUAL 89.
        [_member("DeBreak", "DB16961", "6", "."), _member(DELLY_METHOD, "INV00050937", "2", "89"),
         _member("SVIM-v2.0.0", "svim", "8", "10")],
        [_member(DELLY_METHOD, "INV00050937", "2", "89"), _member("DeBreak", "DB16961", "6", ".")],
    ):
        named = select_representative_sv(group).sv_id
        for member in group:
            if member.sample["SC"] == DELLY_METHOD:
                member.sample = {"SC": "."}
        assert named == select_representative_sv(group).sv_id == "DB16961"


def test_filled_name_is_the_only_difference(tmp_path):
    # Without ##source the record is corrected exactly as if the VCF had
    # declared its own SVMETHOD as ##source: same records, same SC.
    (tmp_path / "a").mkdir()
    (tmp_path / "b").mkdir()
    filled = _correct(tmp_path / "a", _delly_vcf(tmp_path / "a" / "delly.vcf", source_line=None))
    declared = _correct(tmp_path / "b", _delly_vcf(tmp_path / "b" / "delly.vcf",
                                                     source_line=f"##source={DELLY_METHOD}"))

    def data(path):
        return [line for line in path.read_text().splitlines() if line and not line.startswith("#")]

    assert data(filled) == data(declared)
