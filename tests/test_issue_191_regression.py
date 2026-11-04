"""Issue #191: --tra-strand-consistency must be switchable off from the CLI.

The default (orientation consistency on) must stay exactly as before; the new
--no-tra-strand-consistency flag only changes TRA matching.
"""

from pathlib import Path

from typer.testing import CliRunner

from octopusv.cli.cli import app


runner = CliRunner()


def _caller_svcf(path: Path, record_id: str, caller: str, alt: str) -> Path:
    info = ";".join(
        [
            "SVTYPE=TRA",
            "END=5000",
            "SVLEN=.",
            "CHR2=chr2",
            "SUPPORT=8",
            "SVMETHOD=OctopuSV",
            "RTID=.",
            "AF=.",
            "STRAND=.",
            "RNAMES=.",
        ]
    )
    evidence = ":".join(
        ["0/1", ".,.", ".", ".", "60", "TRA", record_id, caller, "N", alt,
         "chr1_1000-chr2_5000"]
    )
    path.write_text(
        "\n".join(
            [
                "##fileformat=VCFv4.2",
                "##SVCFVersion=1.1",
                "##OctopuSV_mode=caller",
                "##contig=<ID=chr1,length=248956422>",
                "##contig=<ID=chr2,length=242193529>",
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE",
                "\t".join(
                    [
                        "chr1", "1000", record_id, "N", alt, "60", "PASS", info,
                        "GT:AD:LN:ST:QV:TY:ID:SC:REF:ALT:CO", evidence,
                    ]
                ),
            ]
        )
        + "\n"
    )
    return path


def _inputs(tmp_path: Path) -> list[str]:
    # Same coordinates, incompatible breakend adjacencies:
    # t[p[ joins chr1-left to chr2-right; [p[t joins chr1-right to chr2-right.
    a = _caller_svcf(tmp_path / "callerA.svcf", "a1", "callerA", "N[chr2:5000[")
    b = _caller_svcf(tmp_path / "callerB.svcf", "b1", "callerB", "[chr2:5000[N")
    return [str(a), str(b)]


def _records(path: Path) -> list[str]:
    return [line for line in path.read_text().splitlines() if not line.startswith("#")]


def _body_without_file_date(path: Path) -> list[str]:
    return [
        line for line in path.read_text().splitlines()
        if not line.startswith("##fileDate=")
    ]


def _merge(inputs: list[str], output: Path, *extra: str):
    args = ["merge", "-i", inputs[0], "-i", inputs[1], "--union", "-o", str(output), *extra]
    return runner.invoke(app, args)


def test_default_keeps_orientation_consistency(tmp_path: Path):
    inputs = _inputs(tmp_path)
    out = tmp_path / "default.svcf"
    result = _merge(inputs, out)
    assert result.exit_code == 0, result.output
    assert len(_records(out)) == 2


def test_explicit_flag_is_identical_to_default(tmp_path: Path):
    inputs = _inputs(tmp_path)
    default_out = tmp_path / "default.svcf"
    explicit_out = tmp_path / "explicit.svcf"
    assert _merge(inputs, default_out).exit_code == 0
    assert _merge(inputs, explicit_out, "--tra-strand-consistency").exit_code == 0
    assert _body_without_file_date(default_out) == _body_without_file_date(explicit_out)


def test_no_flag_disables_orientation_consistency(tmp_path: Path):
    inputs = _inputs(tmp_path)
    out = tmp_path / "relaxed.svcf"
    result = _merge(inputs, out, "--no-tra-strand-consistency")
    assert result.exit_code == 0, result.output
    records = _records(out)
    assert len(records) == 1
    assert "SOURCES=callerA,callerB" in records[0]
