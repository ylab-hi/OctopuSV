from __future__ import annotations

import json
from pathlib import Path

import pytest
from typer.testing import CliRunner

pytest.importorskip("natsort")

from octopusv.cli.cli import app


CALLER_FORMAT = "GT:AD:LN:ST:QV:TY:ID:SC:REF:ALT:CO"


def _invoke(args: list[str]):
    return CliRunner().invoke(app, args)


def _caller_svcf(*, pos: int = 100, record_id: str = "event1", precise: bool = False) -> str:
    end = pos + 100
    precise_token = ";PRECISE" if precise else ""
    block = (
        f"0/1:5,7:100:+-:60:DEL:{record_id}:caller:N:<DEL>:"
        f"chr1_{pos}-chr1_{end}"
    )
    return (
        "##fileformat=VCFv4.2\n"
        "##SVCFVersion=1.1\n"
        "##OctopuSV_mode=caller\n"
        "##contig=<ID=chr1,length=1000000>\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n"
        f"chr1\t{pos}\t{record_id}\tN\t<DEL>\t60\tPASS\t"
        f"SVTYPE=DEL;END={end};SVLEN=100;CHR2=chr1;SUPPORT=5;"
        f"SVMETHOD=OctopuSV;RTID=.;AF=.;STRAND=+-;RNAMES=.{precise_token}\t"
        f"{CALLER_FORMAT}\t{block}\n"
    )


def test_all_registered_commands_render_help():
    commands = [
        "header",
        "validate-svcf",
        "stat",
        "inspect",
        "filter",
        "query",
        "subset",
        "clean",
        "correct",
        "merge",
        "normalize-contigs",
        "svcf2vcf",
        "svcf2bed",
        "svcf2bedpe",
        "benchmark",
        "plot",
        "plot-circos",
        "somatic",
    ]

    top = _invoke(["--help"])
    assert top.exit_code == 0, top.output

    for command in commands:
        result = _invoke([command, "--help"])
        assert result.exit_code == 0, f"{command}: {result.output}"


def test_somatic_wrapper_calls_merge_with_real_none_defaults(tmp_path):
    tumor = tmp_path / "tumor.svcf"
    normal = tmp_path / "normal.svcf"
    output = tmp_path / "somatic.svcf"
    tumor.write_text(_caller_svcf(pos=100, record_id="tumor.event"))
    normal.write_text(_caller_svcf(pos=1000, record_id="normal.event"))

    result = _invoke(
        [
            "somatic",
            "--tumor",
            str(tumor),
            "--normal",
            str(normal),
            "--output-file",
            str(output),
        ]
    )

    assert result.exit_code == 0, result.output
    text = output.read_text()
    assert "##SVCFVersion=1.1\n" in text
    assert "##OctopuSV_mode=multi\n" in text
    assert "tumor.event" in text
    assert "normal.event" not in text


def test_benchmark_preserves_bare_info_flags(tmp_path):
    truth = tmp_path / "truth.svcf"
    call = tmp_path / "call.svcf"
    out_dir = tmp_path / "bench"
    text = _caller_svcf(precise=True)
    truth.write_text(text)
    call.write_text(text)

    result = _invoke(["benchmark", str(truth), str(call), "-o", str(out_dir)])
    assert result.exit_code == 0, result.output

    summary = json.loads((out_dir / "summary.json").read_text())
    assert summary["TP"] == 1
    assert summary["FP"] == 0
    assert summary["FN"] == 0

    for name in ("tp-base.vcf", "tp-call.vcf"):
        record = [
            line
            for line in (out_dir / name).read_text().splitlines()
            if line and not line.startswith("#")
        ][0]
        info_tokens = record.split("\t")[7].split(";")
        assert "PRECISE" in info_tokens
        assert "PRECISE=True" not in info_tokens


def test_clean_missing_external_tools_is_actionable(tmp_path, monkeypatch):
    import octopusv.cli.clean as clean_module

    input_vcf = tmp_path / "input.vcf"
    input_vcf.write_text(
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
    )

    monkeypatch.setattr(clean_module.shutil, "which", lambda _tool: None)
    result = _invoke(["clean", "-i", str(input_vcf), "-o", str(tmp_path / "out.vcf.gz")])

    assert result.exit_code == 1
    assert "required external tool(s) not found" in result.output
    assert "bcftools" in result.output
    assert "bgzip" in result.output
    assert "tabix" in result.output


def test_plot_cli_wires_all_three_standard_plots(tmp_path, monkeypatch):
    import octopusv.cli.plot as plot_module

    stats = tmp_path / "stat.json"
    stats.write_text("{}")
    calls: list[tuple[str, str, bool]] = []

    class FakePlotter:
        def __init__(self, _input):
            pass

        def plot(self, output_prefix, *, save_svg=True):
            calls.append((self.__class__.__name__, str(output_prefix), save_svg))

    class FakeChromosomePlotter(FakePlotter):
        pass

    class FakeTypePlotter(FakePlotter):
        pass

    class FakeSizePlotter(FakePlotter):
        pass

    monkeypatch.setattr(plot_module, "ChromosomePlotter", FakeChromosomePlotter)
    monkeypatch.setattr(plot_module, "TypePlotter", FakeTypePlotter)
    monkeypatch.setattr(plot_module, "SizePlotter", FakeSizePlotter)

    prefix = tmp_path / "plots"
    result = _invoke(["plot", "-i", str(stats), "-o", str(prefix), "--no-svg"])

    assert result.exit_code == 0, result.output
    assert calls == [
        ("FakeChromosomePlotter", f"{prefix}_chromosome_distribution", False),
        ("FakeTypePlotter", f"{prefix}_sv_types", False),
        ("FakeSizePlotter", f"{prefix}_sv_sizes", False),
    ]


def test_plot_circos_missing_optional_dependency_is_actionable(tmp_path, monkeypatch):
    import octopusv.cli.plot_circos as circos_module

    input_svcf = tmp_path / "input.svcf"
    input_svcf.write_text(_caller_svcf())

    real_import_module = circos_module.importlib.import_module

    def fake_import_module(name, *args, **kwargs):
        if name == "pycirclize":
            exc = ModuleNotFoundError("No module named 'pycirclize'")
            exc.name = "pycirclize"
            raise exc
        return real_import_module(name, *args, **kwargs)

    monkeypatch.setattr(circos_module.importlib, "import_module", fake_import_module)

    result = _invoke(
        [
            "plot-circos",
            "-i",
            str(input_svcf),
            "-o",
            str(tmp_path / "circos.png"),
        ]
    )

    assert result.exit_code == 1
    assert "pycirclize is required" in result.output
    assert "pip install pycirclize" in result.output
