from __future__ import annotations

import pytest
from typer.testing import CliRunner

pytest.importorskip("natsort")

from octopusv.cli.cli import app


runner = CliRunner()


def test_issue_179_rejects_intersect_with_min_support(tmp_path):
    """Conflicting merge strategies must fail instead of silently ignoring one."""
    input_a = tmp_path / "a.svcf"
    input_b = tmp_path / "b.svcf"
    output = tmp_path / "merged.svcf"

    # The strategy guard should run before input parsing, so these files do not
    # need to contain valid SVCF records.
    input_a.write_text("", encoding="utf-8")
    input_b.write_text("", encoding="utf-8")

    result = runner.invoke(
        app,
        [
            "merge",
            "-i",
            str(input_a),
            "-i",
            str(input_b),
            "-o",
            str(output),
            "--intersect",
            "--min-support",
            "3",
        ],
    )

    assert result.exit_code == 1
    assert "Conflicting merge strategies" in result.output
    assert "--intersect" in result.output
    assert "--min-support/--max-support" in result.output
    assert not output.exists()
