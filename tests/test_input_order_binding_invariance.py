from __future__ import annotations

from itertools import permutations
from pathlib import Path

import pytest
from typer.testing import CliRunner

pytest.importorskip("natsort")

from octopusv.cli.cli import app
from octopusv.utils.svcf_sample_parser import parse_svcf_sample_block


runner = CliRunner()
CALLER_FORMAT = "GT:AD:LN:ST:QV:TY:ID:SC:REF:ALT:CO"


def _write_input(path: Path, label: str, record_id: str, gt: str) -> None:
    block = (
        f"{gt}:5,5:100:.:60:DEL:{record_id}:caller{label}:N:<DEL>:"
        "chr1_100-chr1_200"
    )
    path.write_text(
        "##fileformat=VCFv4.2\n"
        "##SVCFVersion=1.1\n"
        "##OctopuSV_mode=caller\n"
        "##contig=<ID=chr1,length=1000000>\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n"
        f"chr1\t100\t{record_id}\tN\t<DEL>\t60\tPASS\t"
        "SVTYPE=DEL;END=200;SVLEN=100;CHR2=chr1;SUPPORT=5;"
        "SVMETHOD=OctopuSV;RTID=.;AF=.;STRAND=.;RNAMES=.\t"
        f"{CALLER_FORMAT}\t{block}\n"
    )


def _invoke(args: list[str]) -> None:
    result = runner.invoke(app, args)
    assert result.exit_code == 0, result.output


def _one_record(path: Path) -> list[str]:
    records = [
        line.split("\t")
        for line in path.read_text().splitlines()
        if line and not line.startswith("#")
    ]
    assert len(records) == 1
    return records[0]


def _parse_info(text: str) -> dict[str, str | bool]:
    out: dict[str, str | bool] = {}
    for item in text.split(";"):
        if "=" in item:
            key, value = item.split("=", 1)
            out[key] = value
        elif item:
            out[item] = True
    return out


def _caller_binding_signature(path: Path):
    row = _one_record(path)
    info = _parse_info(row[7])
    sources = str(info["SOURCES"]).split(",")
    source_ids = str(info["SOURCE_IDS"]).split(",")
    blocks = [parse_svcf_sample_block(row[8], block) for block in row[9:]]
    assert len(sources) == len(source_ids) == len(blocks)
    triples = {
        (source, source_id, str(block["ID"]), str(block["SC"]), str(block["GT"]))
        for source, source_id, block in zip(sources, source_ids, blocks)
    }
    invariant = (row[0], row[1], row[3], row[4], info["SVTYPE"], info["END"], info["SVLEN"])
    return invariant, triples


def _sample_binding_signature(path: Path):
    lines = path.read_text().splitlines()
    chrom = next(line for line in lines if line.startswith("#CHROM"))
    labels = chrom.split("\t")[9:]
    row = _one_record(path)
    assert len(labels) == len(row[9:])
    mapping = {label: block for label, block in zip(labels, row[9:])}
    info = _parse_info(row[7])
    invariant = (row[0], row[1], row[3], row[4], info["SVTYPE"], info["END"], info["SVLEN"])
    return invariant, mapping


def _vcf_sample_gt_signature(path: Path):
    lines = path.read_text().splitlines()
    chrom = next(line for line in lines if line.startswith("#CHROM"))
    labels = chrom.split("\t")[9:]
    row = _one_record(path)
    fmt = row[8].split(":")
    gt_index = fmt.index("GT")
    mapping = {
        label: block.split(":")[gt_index]
        for label, block in zip(labels, row[9:])
    }
    return (row[0], row[1], row[3], row[4]), mapping


def test_three_input_permutations_preserve_binding_relationships(tmp_path):
    inputs: dict[str, Path] = {}
    for label, record_id, gt in [("A", "idA", "0/1"), ("B", "idB", "1/1"), ("C", "idC", "1/.")]:
        path = tmp_path / f"{label}.svcf"
        _write_input(path, label, record_id, gt)
        inputs[label] = path

    orders = [("A", "B", "C"), ("B", "C", "A"), ("C", "A", "B")]
    caller_signatures = []
    sample_signatures = []
    vcf_signatures = []
    subset_signatures = []

    for index, order in enumerate(orders):
        input_args: list[str] = []
        for label in order:
            input_args += ["-i", str(inputs[label])]

        caller_out = tmp_path / f"caller_{index}.svcf"
        _invoke([
            "merge", *input_args, "-o", str(caller_out), "--mode", "caller",
            "--caller-names", ",".join(order), "--union",
        ])
        caller_signatures.append(_caller_binding_signature(caller_out))

        sample_out = tmp_path / f"sample_{index}.svcf"
        _invoke([
            "merge", *input_args, "-o", str(sample_out), "--mode", "sample",
            "--sample-names", ",".join(order), "--union",
        ])
        sample_signatures.append(_sample_binding_signature(sample_out))

        vcf_out = tmp_path / f"sample_{index}.vcf"
        _invoke(["svcf2vcf", "-i", str(sample_out), "-o", str(vcf_out)])
        vcf_signatures.append(_vcf_sample_gt_signature(vcf_out))

        subset_out = tmp_path / f"subset_{index}.svcf"
        _invoke(["subset", "-i", str(sample_out), "-o", str(subset_out), "--sample", "A"])
        subset_signatures.append(_sample_binding_signature(subset_out))

    # Record-level geometry must be invariant to input order.
    assert len({sig[0] for sig in caller_signatures}) == 1
    assert len({sig[0] for sig in sample_signatures}) == 1

    # Caller mode may order SOURCES by input order, but the source/id/evidence
    # relationship itself must remain identical.
    assert caller_signatures[0][1] == caller_signatures[1][1] == caller_signatures[2][1]

    # Sample columns may follow input order, but label -> data binding must not.
    assert sample_signatures[0][1] == sample_signatures[1][1] == sample_signatures[2][1]
    assert vcf_signatures[0] == vcf_signatures[1] == vcf_signatures[2]
    assert subset_signatures[0] == subset_signatures[1] == subset_signatures[2]
