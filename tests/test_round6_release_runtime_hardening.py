import builtins
import json
import os
import re
import subprocess
import sys
from pathlib import Path
from types import SimpleNamespace

import pytest

from octopusv import __version__
from octopusv.bencher.sv_bencher import SVBencher


REPO_ROOT = Path(__file__).resolve().parents[1]


def _subprocess_env() -> dict[str, str]:
    env = os.environ.copy()
    src = REPO_ROOT / "src"
    package_root = src if src.exists() else REPO_ROOT
    current = env.get("PYTHONPATH", "")
    env["PYTHONPATH"] = str(package_root) + (os.pathsep + current if current else "")
    return env


def _run_python(*args: str) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        [sys.executable, *args],
        cwd=REPO_ROOT,
        env=_subprocess_env(),
        capture_output=True,
        text=True,
        check=False,
    )


def test_cli_version_outputs_only_package_version():
    result = _run_python("-m", "octopusv", "--version")

    assert result.returncode == 0, result.stderr
    assert result.stdout.strip() == __version__
    assert result.stderr == ""


def test_cli_without_command_preserves_missing_command_behavior():
    result = _run_python("-m", "octopusv")

    assert result.returncode == 2
    assert "Missing command" in result.stderr


def test_pyproject_version_matches_package_version_without_tomllib():
    text = (REPO_ROOT / "pyproject.toml").read_text(encoding="utf-8")
    tool_poetry = text.split("[tool.poetry]", 1)[1].split("[", 1)[0]
    match = re.search(r'^version\s*=\s*["\']([^"\']+)["\']\s*$', tool_poetry, re.MULTILINE)

    assert match is not None
    assert match.group(1) == __version__



def test_library_imports_do_not_configure_root_logging():
    script = (
        "import logging\n"
        "import octopusv.converter.base\n"
        "import octopusv.formatter.svcf_to_vcf_converter\n"
        "import octopusv.report.generator\n"
        "import octopusv.cli.bench\n"
        "raise SystemExit(1 if logging.getLogger().handlers else 0)\n"
    )
    result = _run_python("-c", script)

    assert result.returncode == 0, result.stderr

def test_cli_logging_uses_stderr_and_does_not_pollute_json_stdout():
    script = (
        "import json, logging\n"
        "import octopusv.cli.cli\n"
        "logging.getLogger('octopusv.round6').warning('round6-warning')\n"
        "print(json.dumps({'ok': True}))\n"
    )
    result = _run_python("-c", script)

    assert result.returncode == 0, result.stderr
    assert json.loads(result.stdout) == {"ok": True}
    assert "round6-warning" not in result.stdout
    assert "round6-warning" in result.stderr


def test_importing_cli_does_not_import_matplotlib():
    script = (
        "import sys\n"
        "import octopusv.cli.cli\n"
        "raise SystemExit(1 if 'matplotlib' in sys.modules else 0)\n"
    )
    result = _run_python("-c", script)

    assert result.returncode == 0, result.stderr



def test_plot_backends_load_matplotlib_only_when_requested():
    script = (
        "import sys\n"
        "import octopusv.cli.cli\n"
        "assert 'matplotlib' not in sys.modules\n"
        "import octopusv.cli.plot as plot_module\n"
        "plot_module._load_plotters()\n"
        "raise SystemExit(0 if 'matplotlib' in sys.modules else 1)\n"
    )
    result = _run_python("-c", script)

    assert result.returncode == 0, result.stderr


def test_upset_plotter_still_generates_output_when_explicitly_loaded(tmp_path):
    from octopusv.merger.upset_plotter import UpSetPlotter

    event = SimpleNamespace(source_file="a.svcf")
    output = tmp_path / "upset.png"
    UpSetPlotter([event], ["a.svcf"]).plot(str(output))

    assert output.exists()
    assert output.stat().st_size > 0

def test_sequence_comparison_missing_levenshtein_fails_before_benchmark(monkeypatch, tmp_path):
    real_import = builtins.__import__

    def guarded_import(name, *args, **kwargs):
        if name == "Levenshtein":
            raise ImportError("simulated missing Levenshtein")
        return real_import(name, *args, **kwargs)

    monkeypatch.setattr(builtins, "__import__", guarded_import)

    with pytest.raises(RuntimeError, match="requires the Levenshtein package"):
        SVBencher(
            tmp_path / "truth.svcf",
            tmp_path / "call.svcf",
            tmp_path / "out",
            enable_sequence_comparison=True,
        )


def test_sequence_comparison_uses_loaded_ratio_without_fallback(monkeypatch, tmp_path):
    calls: list[tuple[str, str]] = []

    def fake_ratio(left: str, right: str) -> float:
        calls.append((left, right))
        return 0.625

    monkeypatch.setattr(SVBencher, "_load_sequence_ratio", staticmethod(lambda: fake_ratio))
    bencher = SVBencher(
        tmp_path / "truth.svcf",
        tmp_path / "call.svcf",
        tmp_path / "out",
        enable_sequence_comparison=True,
    )

    truth = SimpleNamespace(alt_seq="AACCGG")
    call = SimpleNamespace(alt_seq="AACCTT")

    assert bencher._calculate_sequence_similarity(truth, call) == 0.625
    assert calls == [("AACCGG", "AACCTT")]


def test_pyproject_runtime_dependency_cleanup():
    text = (REPO_ROOT / "pyproject.toml").read_text(encoding="utf-8")
    runtime = text.split("[tool.poetry.dependencies]", 1)[1].split("[", 1)[0]
    dev = text.split("[tool.poetry.group.dev.dependencies]", 1)[1].split("[", 1)[0]

    assert re.search(r"^numpy\s*=", runtime, re.MULTILINE)
    assert re.search(r"^levenshtein\s*=", runtime, re.MULTILINE)
    assert not re.search(r"^pytest-cov\s*=", runtime, re.MULTILINE)
    assert not re.search(r"^seaborn\s*=", runtime, re.MULTILINE)
    assert not re.search(r"^loguru\s*=", runtime, re.MULTILINE)
    # Dependency ranges were finalized after clean-environment validation:
    # minimum: Typer 0.12.4 / Click 8.0.0 / Rich 13.7.1;
    # latest tested: Typer 0.27.2 / Click 8.5.0 / Rich 15.0.0.
    assert re.search(r'^typer\s*=\s*["\']>=0\.12\.4,<0\.28["\']', runtime, re.MULTILINE)
    assert re.search(r'^rich\s*=\s*["\']>=13\.7\.1,<16["\']', runtime, re.MULTILINE)
    assert not re.search(r"^click\s*=", runtime, re.MULTILINE)
    assert re.search(r"^pytest-cov\s*=", dev, re.MULTILINE)


def test_pyproject_has_project_urls():
    text = (REPO_ROOT / "pyproject.toml").read_text(encoding="utf-8")
    tool_poetry = text.split("[tool.poetry]", 1)[1].split("[tool.poetry.dependencies]", 1)[0]

    assert 'homepage = "https://github.com/ylab-hi/OctopuSV"' in tool_poetry
    assert 'repository = "https://github.com/ylab-hi/OctopuSV"' in tool_poetry
    assert 'documentation = "https://github.com/ylab-hi/OctopuSV#readme"' in tool_poetry
