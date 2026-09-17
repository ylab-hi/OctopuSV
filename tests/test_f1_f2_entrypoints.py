from __future__ import annotations

import gzip
import json
from pathlib import Path

import pytest
from typer.testing import CliRunner

pytest.importorskip("natsort")

from octopusv.cli.cli import app


CALLER_FORMAT = "GT:AD:LN:ST:QV:TY:ID:SC:REF:ALT:CO"


def _base_svcf() -> str:
    block = "0/1:5,7:50:+-:60:DEL:caller.1:caller:N:<DEL>:chr1_100-chr1_150"
    return (
        "##fileformat=VCFv4.2\n"
        "##SVCFVersion=1.1\n"
        "##OctopuSV_mode=caller\n"
        "##contig=<ID=chr1,length=1000>\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n"
        "chr1\t100\tevent1\tN\t<DEL>\t60\tPASS\t"
        "SVTYPE=DEL;END=150;SVLEN=50;CHR2=chr1;SUPPORT=5;SVMETHOD=OctopuSV;"
        "RTID=.;AF=.;STRAND=+-;RNAMES=.;SOURCES=caller;SOURCE_IDS=caller.1\t"
        f"{CALLER_FORMAT}\t{block}\n"
    )



def _raw_vcf() -> str:
    return (
        "##fileformat=VCFv4.2\n"
        "##source=TestCaller\n"
        "##contig=<ID=chr1,length=1000>\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n"
        "chr1\t100\tr1\tN\t<DEL>\t60\tPASS\t"
        "SVTYPE=DEL;END=150;SVLEN=-50\tGT:AD\t0/1:5,7\n"
    )

def _write_variant(path: Path, kind: str, text: str) -> None:
    if kind in {"gzip_suffix", "gzip_no_suffix"}:
        with gzip.open(path, "wt", encoding="utf-8") as handle:
            handle.write(text)
    elif kind == "bom":
        path.write_bytes(b"\xef\xbb\xbf" + text.encode("utf-8"))
    else:
        path.write_text(text)


def _normalize_text(text: str) -> str:
    return "\n".join(
        line for line in text.splitlines() if not line.startswith("##fileDate=")
    ) + ("\n" if text.endswith("\n") else "")


def _invoke(args: list[str]):
    return CliRunner().invoke(app, args)


READ_VARIANTS = ["plain", "gzip_suffix", "gzip_no_suffix", "bom"]


@pytest.mark.parametrize("kind", READ_VARIANTS)
@pytest.mark.parametrize(
    "command,extra_args",
    [
        ("validate-svcf", ["--json"]),
        ("header", ["--json"]),
        ("inspect", ["--id", "event1", "--json"]),
        ("stat", ["--json"]),
    ],
)
def test_read_only_entrypoints_accept_compression_and_bom(tmp_path, kind, command, extra_args):
    suffix = ".gz" if kind == "gzip_suffix" else ".svcf"
    path = tmp_path / f"input_{kind}{suffix}"
    _write_variant(path, kind, _base_svcf())

    result = _invoke([command, "-i", str(path), *extra_args])
    assert result.exit_code == 0, result.output

    # JSON-producing commands should still return structurally valid output.
    payload = json.loads(result.stdout)
    if isinstance(payload, dict):
        payload.pop("input", None)
    assert payload is not None


@pytest.mark.parametrize("kind", READ_VARIANTS)
def test_correct_is_input_representation_invariant(tmp_path, kind):
    plain = tmp_path / "plain.vcf"
    plain.write_text(_raw_vcf())
    variant_suffix = ".gz" if kind == "gzip_suffix" else ".vcf"
    variant = tmp_path / f"variant_{kind}{variant_suffix}"
    _write_variant(variant, kind, _raw_vcf())

    baseline_out = tmp_path / "baseline.svcf"
    observed_out = tmp_path / f"observed_{kind}.svcf"
    r1 = _invoke(["correct", "-i", str(plain), "-o", str(baseline_out)])
    r2 = _invoke(["correct", "-i", str(variant), "-o", str(observed_out)])
    assert r1.exit_code == 0, r1.output
    assert r2.exit_code == 0, r2.output
    assert _normalize_text(observed_out.read_text()) == _normalize_text(baseline_out.read_text())


@pytest.mark.parametrize("kind", READ_VARIANTS)
@pytest.mark.parametrize(
    "command,extra_args,extension",
    [
        ("filter", ["--svtype", "DEL"], ".svcf"),
        ("query", ["--region", "chr1:1-500"], ".svcf"),
        ("subset", ["--mode", "caller", "--caller", "caller"], ".svcf"),
        ("normalize-contigs", ["--style", "chr"], ".svcf"),
        ("svcf2vcf", [], ".vcf"),
        ("svcf2bed", [], ".bed"),
        ("svcf2bedpe", [], ".bedpe"),
    ],
)
def test_file_writing_entrypoints_are_input_representation_invariant(
    tmp_path, kind, command, extra_args, extension
):
    plain = tmp_path / "plain.svcf"
    plain.write_text(_base_svcf())
    variant_suffix = ".gz" if kind == "gzip_suffix" else ".svcf"
    variant = tmp_path / f"variant_{kind}{variant_suffix}"
    _write_variant(variant, kind, _base_svcf())

    baseline_out = tmp_path / f"baseline_{command}{extension}"
    observed_out = tmp_path / f"observed_{command}_{kind}{extension}"

    r1 = _invoke([command, "-i", str(plain), "-o", str(baseline_out), *extra_args])
    r2 = _invoke([command, "-i", str(variant), "-o", str(observed_out), *extra_args])
    assert r1.exit_code == 0, r1.output
    assert r2.exit_code == 0, r2.output
    assert _normalize_text(observed_out.read_text()) == _normalize_text(baseline_out.read_text())


def _truncated_svcf() -> str:
    return _base_svcf() + "chr1\t300\ttruncated\n"


@pytest.mark.parametrize(
    "command,extra_args,needs_output",
    [
        ("validate-svcf", [], False),
        ("inspect", ["--id", "event1"], False),
        ("stat", ["--json"], False),
        ("filter", ["--svtype", "DEL"], True),
        ("query", ["--region", "chr1:1-500"], True),
        ("subset", ["--mode", "caller", "--caller", "caller"], True),
        ("normalize-contigs", ["--style", "chr"], True),
        ("svcf2vcf", [], True),
        ("svcf2bed", [], True),
        ("svcf2bedpe", [], True),
    ],
)
def test_record_consumers_reject_truncated_rows_with_line_number(
    tmp_path, command, extra_args, needs_output
):
    input_path = tmp_path / "truncated.svcf"
    input_path.write_text(_truncated_svcf())
    args = [command, "-i", str(input_path)]
    output = tmp_path / f"{command}.out"
    if needs_output:
        output.write_text("ORIGINAL\n")
        args += ["-o", str(output)]
    args += extra_args

    result = _invoke(args)
    assert result.exit_code != 0
    combined = result.output or ""
    if result.exception is not None:
        combined += " " + str(result.exception)
    assert "line 7" in combined

    # Commands converted to atomic streaming output must preserve an existing
    # destination. BED/BEDPE materialize conversion before opening their output,
    # so they also leave this sentinel untouched when parsing fails.
    if needs_output:
        assert output.read_text() == "ORIGINAL\n"
