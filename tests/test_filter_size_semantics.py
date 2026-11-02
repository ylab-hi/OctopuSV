"""octopusv filter --min-size: records without a defined length are kept.

Translocations / breakends have no linear size, and an insertion without
SVLEN has an unknown length.  Neither is "smaller than" a threshold
(issue #193, same rule as correct --min-svlen).
"""
from __future__ import annotations

import json

import pytest
from typer.testing import CliRunner

pytest.importorskip("natsort")

from octopusv.cli.cli import app

runner = CliRunner()


def _svcf(tmp_path):
    records = [
        "chr1\t100\ttra\tN\t<TRA>\t.\tPASS\tSVTYPE=TRA;CHR2=chr2;END=500;SVLEN=0;SUPPORT=8\tGT\t0/1",
        "chr1\t300\tbnd\tN\tN[chr2:900[\t.\tPASS\tSVTYPE=BND;SUPPORT=6\tGT\t0/1",
        "chr1\t1000\tdel.small\tN\t<DEL>\t.\tPASS\tSVTYPE=DEL;END=1030;SVLEN=-30;SUPPORT=8\tGT\t0/1",
        "chr1\t2000\tdel.big\tN\t<DEL>\t.\tPASS\tSVTYPE=DEL;END=2100;SVLEN=-100;SUPPORT=8\tGT\t0/1",
        "chr1\t3000\tins.small\tN\t<INS>\t.\tPASS\tSVTYPE=INS;END=3000;SVLEN=20;SUPPORT=8\tGT\t0/1",
        "chr1\t4000\tins.unknown\tN\t<INS>\t.\tPASS\tSVTYPE=INS;END=4020;SUPPORT=8\tGT\t0/1",
    ]
    vcf = tmp_path / "in.vcf"
    vcf.write_text(
        "##fileformat=VCFv4.2\n##source=DeBreak\n"
        "##contig=<ID=chr1,length=1000000>\n##contig=<ID=chr2,length=1000000>\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n"
        + "\n".join(records) + "\n")
    svcf = tmp_path / "in.svcf"
    result = runner.invoke(app, ["correct", "-i", str(vcf), "-o", str(svcf)])
    assert result.exit_code == 0, result.output
    return svcf


def _ids(path):
    return {line.split("\t")[2] for line in path.read_text().splitlines() if line and not line.startswith("#")}


def test_min_size_keeps_records_without_defined_length(tmp_path):
    svcf = _svcf(tmp_path)
    out = tmp_path / "out.svcf"
    result = runner.invoke(app, ["filter", "-i", str(svcf), "-o", str(out), "--min-size", "50",
                                 "--json-summary"])
    assert result.exit_code == 0, result.output
    assert _ids(out) == {"tra", "bnd", "del.big", "ins.unknown"}
    data = json.loads(result.stdout[result.stdout.index("{"):])
    assert data["size_filter_not_applied"] == {"TRA_or_BND_no_linear_size": 2, "INS_length_unknown": 1}
    assert data["excluded_by_reason"] == {"size_below_min": 2}


def test_exclude_svtype_still_removes_translocations_on_request(tmp_path):
    svcf = _svcf(tmp_path)
    out = tmp_path / "out.svcf"
    result = runner.invoke(app, ["filter", "-i", str(svcf), "-o", str(out), "--min-size", "50",
                                 "--exclude-svtype", "TRA,BND"])
    assert result.exit_code == 0, result.output
    assert _ids(out) == {"del.big", "ins.unknown"}
