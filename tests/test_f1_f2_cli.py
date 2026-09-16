from __future__ import annotations

import gzip
from pathlib import Path

import pytest
from typer.testing import CliRunner

pytest.importorskip("natsort")

from octopusv.cli.cli import app


CALLER_FORMAT = "GT:AD:LN:ST:QV:TY:ID:SC:REF:ALT:CO"


def _write_gzip(path: Path, text: str) -> None:
    with gzip.open(path, "wt", encoding="utf-8") as handle:
        handle.write(text)


def _record(pos: int, record_id: str, source: str, source_id: str | None = None) -> str:
    end = pos + 50
    source_id = source_id or record_id
    block = (
        f"0/1:5,7:50:+-:60:DEL:{source_id}:{source}:N:<DEL>:"
        f"chr1_{pos}-chr1_{end}"
    )
    return (
        f"chr1\t{pos}\t{record_id}\tN\t<DEL>\t60\tPASS\t"
        f"SVTYPE=DEL;END={end};SVLEN=50;CHR2=chr1;SUPPORT=5;"
        f"SVMETHOD=OctopuSV;RTID=.;AF=.;STRAND=+-;RNAMES=.\t"
        f"{CALLER_FORMAT}\t{block}\n"
    )


def _caller_svcf_text(records: list[str]) -> str:
    return (
        "##fileformat=VCFv4.2\n"
        "##SVCFVersion=1.1\n"
        "##OctopuSV_mode=caller\n"
        "##contig=<ID=chr1,length=10000000>\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n"
        + "".join(records)
    )


def _normalize_svcf(text: str) -> str:
    return "\n".join(
        line for line in text.splitlines() if not line.startswith("##fileDate=")
    ) + "\n"


def _invoke(args: list[str]):
    return CliRunner().invoke(app, args)


def _make_two_inputs(directory: Path, *, gzip_data: bool = False) -> tuple[Path, Path]:
    directory.mkdir(parents=True, exist_ok=True)
    a = directory / "a.svcf"
    b = directory / "b.svcf"
    a_text = _caller_svcf_text([_record(100, "a.1", "callerA")])
    b_text = _caller_svcf_text([_record(100, "b.1", "callerB")])
    if gzip_data:
        _write_gzip(a, a_text)
        _write_gzip(b, b_text)
    else:
        a.write_text(a_text)
        b.write_text(b_text)
    return a, b


def _run_merge(tmp_path: Path, name: str, args: list[str]) -> str:
    output = tmp_path / f"{name}.svcf"
    result = _invoke(["merge", *args, "--union", "-o", str(output)])
    assert result.exit_code == 0, result.output
    return _normalize_svcf(output.read_text())


@pytest.mark.parametrize("mode", ["caller", "sample"])
def test_merge_input_representation_invariance(mode, tmp_path):
    plain_a, plain_b = _make_two_inputs(tmp_path / "plain")

    gzip_dir = tmp_path / "gzip"
    gzip_dir.mkdir()
    gz_a = gzip_dir / "a.svcf.gz"
    gz_b = gzip_dir / "b.svcf.gz"
    _write_gzip(gz_a, plain_a.read_text())
    _write_gzip(gz_b, plain_b.read_text())

    magic_a, magic_b = _make_two_inputs(tmp_path / "magic", gzip_data=True)

    links = tmp_path / "links"
    links.mkdir()
    link_a = links / "a.svcf"
    link_b = links / "b.svcf"
    link_a.symlink_to(plain_a)
    link_b.symlink_to(plain_b)
    lists = tmp_path / "lists"
    lists.mkdir()
    input_list = lists / "inputs.tsv"
    input_list.write_text("../links/a.svcf\n../links/b.svcf\n")

    expected = _run_merge(
        tmp_path,
        f"{mode}_positional",
        [str(plain_a), str(plain_b), "--mode", mode],
    )

    variants = [
        ["-i", str(plain_a), "-i", str(plain_b), "--mode", mode],
        ["--input-list", str(input_list), "--mode", mode],
        ["-i", str(gz_a), "-i", str(gz_b), "--mode", mode],
        ["-i", str(magic_a), "-i", str(magic_b), "--mode", mode],
    ]
    for index, args in enumerate(variants):
        observed = _run_merge(tmp_path, f"{mode}_variant_{index}", args)
        assert observed == expected


@pytest.mark.parametrize("mode,label_option", [("caller", "--caller-names"), ("sample", "--sample-names")])
def test_labeled_input_list_matches_explicit_name_vector(mode, label_option, tmp_path):
    a, b = _make_two_inputs(tmp_path / "inputs")
    direct = _run_merge(
        tmp_path,
        f"{mode}_direct_labels",
        ["-i", str(a), "-i", str(b), "--mode", mode, label_option, "A,B"],
    )

    list_path = tmp_path / "labeled.tsv"
    list_path.write_text(f"{a}\tA\n{b}\tB\n")
    listed = _run_merge(
        tmp_path,
        f"{mode}_list_labels",
        ["--input-list", str(list_path), "--mode", mode],
    )
    assert listed == direct


def test_explicit_labels_reject_mixed_input_mechanisms(tmp_path):
    a, b = _make_two_inputs(tmp_path / "inputs")
    result = _invoke(
        [
            "merge",
            str(b),
            f"--input-file={a}",
            "--mode",
            "caller",
            "--caller-names",
            "A,B",
            "--union",
            "-o",
            str(tmp_path / "out.svcf"),
        ]
    )
    assert result.exit_code != 0
    assert "cannot be combined with multiple input mechanisms" in result.output
    assert not (tmp_path / "out.svcf").exists()


def test_labeled_input_list_rejects_sample_names(tmp_path):
    a, b = _make_two_inputs(tmp_path / "inputs")
    list_path = tmp_path / "labeled.tsv"
    list_path.write_text(f"{a}\tA\n{b}\tB\n")
    result = _invoke(
        [
            "merge",
            "--input-list",
            str(list_path),
            "--mode",
            "sample",
            "--sample-names",
            "X,Y",
            "--union",
            "-o",
            str(tmp_path / "out.svcf"),
        ]
    )
    assert result.exit_code != 0
    assert "cannot be combined" in result.output
    assert not (tmp_path / "out.svcf").exists()


def test_correct_cli_accepts_gzip_without_suffix(tmp_path):
    raw = (
        "##fileformat=VCFv4.2\n"
        "##source=TestCaller\n"
        "##contig=<ID=chr1,length=1000>\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n"
        "chr1\t100\tr1\tN\t<DEL>\t60\tPASS\t"
        "SVTYPE=DEL;END=150;SVLEN=-50\tGT:AD\t0/1:5,7\n"
    )
    input_vcf = tmp_path / "input.vcf"
    _write_gzip(input_vcf, raw)
    output = tmp_path / "output.svcf"

    result = _invoke(["correct", "-i", str(input_vcf), "-o", str(output)])
    assert result.exit_code == 0, result.output
    text = output.read_text()
    assert "##SVCFVersion=1.1\n" in text
    assert "##OctopuSV_mode=caller\n" in text


@pytest.mark.parametrize("mode", ["caller", "sample"])
def test_merge_real_writer_failure_is_atomic(mode, tmp_path):
    input_path = tmp_path / "a.svcf"
    input_path.write_text(
        _caller_svcf_text(
            [
                _record(100, "event.good", "callerA", "good.id"),
                _record(1000, "bad id", "callerA", "bad id"),
            ]
        )
    )
    output = tmp_path / "out.svcf"
    output.write_text("ORIGINAL\n")

    result = _invoke(
        [
            "merge",
            "-i",
            str(input_path),
            "--mode",
            mode,
            "--union",
            "-o",
            str(output),
        ]
    )
    assert result.exit_code != 0
    assert "bad id" in result.output
    assert output.read_text() == "ORIGINAL\n"
    assert not list(tmp_path.glob(".out.svcf.*.tmp"))



def test_merge_truncated_input_preserves_existing_output(tmp_path):
    input_path = tmp_path / "truncated.svcf"
    input_path.write_text(
        _caller_svcf_text([_record(100, "event.good", "callerA")])
        + "chr1\t300\ttruncated\n"
    )
    output = tmp_path / "out.svcf"
    output.write_text("ORIGINAL\n")

    result = _invoke(
        [
            "merge",
            "-i",
            str(input_path),
            "--mode",
            "caller",
            "--union",
            "-o",
            str(output),
        ]
    )
    assert result.exit_code != 0
    assert "line 7" in result.output
    assert output.read_text() == "ORIGINAL\n"
    assert not list(tmp_path.glob(".out.svcf.*.tmp"))


@pytest.mark.parametrize("label_source", ["default", "caller_names", "sample_names", "input_list"])
def test_invalid_labels_from_all_user_entrypoints_fail_before_event_parse(
    label_source, tmp_path, monkeypatch
):
    import octopusv.cli.merge as merge_module

    valid = tmp_path / "valid.svcf"
    valid.write_text(_caller_svcf_text([_record(100, "event.good", "callerA")]))
    output = tmp_path / "out.svcf"

    def parse_must_not_run(self):
        raise AssertionError("event parsing must not start before label validation")

    monkeypatch.setattr(merge_module.SVCFFileEventCreator, "parse", parse_must_not_run)

    if label_source == "default":
        invalid_path = tmp_path / "patient 1.svcf"
        invalid_path.write_text(valid.read_text())
        args = ["merge", "-i", str(invalid_path), "--mode", "caller"]
    elif label_source == "caller_names":
        args = [
            "merge", "-i", str(valid), "--mode", "caller",
            "--caller-names", "bad name",
        ]
    elif label_source == "sample_names":
        args = [
            "merge", "-i", str(valid), "--mode", "sample",
            "--sample-names", "bad name",
        ]
    else:
        list_path = tmp_path / "inputs.tsv"
        list_path.write_text(f"{valid}\tbad name\n")
        args = ["merge", "--input-list", str(list_path), "--mode", "sample"]

    result = _invoke([*args, "--union", "-o", str(output)])
    assert result.exit_code != 0
    assert "Invalid" in result.output and "label" in result.output
    assert "patient_1" not in result.output
    assert not output.exists()

def test_bnd_export_contract_cli(tmp_path):
    block = (
        "0/1:5,7:.:+-:60:BND:bnd.1:caller:N:N]chr2:250]:"
        "chr1_100-chr2_250"
    )
    text = (
        "##fileformat=VCFv4.2\n"
        "##SVCFVersion=1.1\n"
        "##OctopuSV_mode=caller\n"
        "##contig=<ID=chr1,length=1000>\n"
        "##contig=<ID=chr2,length=1000>\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n"
        "chr1\t100\tbnd.1\tN\tN]chr2:250]\t60\tPASS\t"
        "SVTYPE=BND;END=250;SVLEN=.;CHR2=chr2;SUPPORT=5;SVMETHOD=OctopuSV;"
        "RTID=.;AF=.;STRAND=+-;RNAMES=.\t"
        f"{CALLER_FORMAT}\t{block}\n"
    )
    input_path = tmp_path / "bnd.svcf"
    input_path.write_text(text)
    bed = tmp_path / "out.bed"
    bedpe = tmp_path / "out.bedpe"

    r1 = _invoke(["svcf2bed", "-i", str(input_path), "-o", str(bed)])
    r2 = _invoke(["svcf2bedpe", "-i", str(input_path), "-o", str(bedpe)])
    assert r1.exit_code == 0, r1.output
    assert r2.exit_code == 0, r2.output

    bed_fields = bed.read_text().splitlines()[-1].split("\t")
    assert bed_fields[:3] == ["chr1", "99", "100"]
    assert "_BND_" in bed_fields[3]

    bedpe_fields = bedpe.read_text().splitlines()[-1].split("\t")
    assert bedpe_fields[:6] == ["chr1", "99", "100", "chr2", "249", "250"]
