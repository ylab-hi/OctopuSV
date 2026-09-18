from pathlib import Path

import pytest
from typer.testing import CliRunner

pytest.importorskip("natsort")

from octopusv.cli.cli import app


runner = CliRunner()

CALLER_FORMAT = "GT:AD:LN:ST:QV:TY:ID:SC:REF:ALT:CO"


def _invoke(args: list[str]):
    result = runner.invoke(app, args)
    assert result.exit_code == 0, (
        f"Command failed: octopusv {' '.join(args)}\n"
        f"output:\n{result.output}\n"
        f"exception: {result.exception!r}"
    )
    return result


def _write_del(
    path: Path,
    *,
    record_id: str,
    pos: int,
    end: int,
    source: str,
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)

    svlen = end - pos

    info = ";".join(
        [
            "SVTYPE=DEL",
            f"END={end}",
            f"SVLEN={svlen}",
            "CHR2=chr1",
            "SUPPORT=10",
            "SVMETHOD=OctopuSV",
            "RTID=.",
            "AF=.",
            "STRAND=.",
            "RNAMES=.",
        ]
    )

    block = ":".join(
        [
            "0/1",
            ".,.",
            str(svlen),
            ".",
            "60",
            "DEL",
            record_id,
            source,
            "N",
            "<DEL>",
            f"chr1_{pos}-chr1_{end}",
        ]
    )

    path.write_text(
        "\n".join(
            [
                "##fileformat=VCFv4.2",
                "##SVCFVersion=1.1",
                "##OctopuSV_mode=caller",
                "##contig=<ID=chr1,length=1000000>",
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1",
                (
                    f"chr1\t{pos}\t{record_id}\tN\t<DEL>\t60\tPASS\t"
                    f"{info}\t{CALLER_FORMAT}\t{block}"
                ),
            ]
        )
        + "\n",
        encoding="utf-8",
    )


def _record_count(path: Path) -> int:
    return sum(
        1
        for line in path.read_text(encoding="utf-8").splitlines()
        if line and not line.startswith("#")
    )


def _merge_pair(
    tmp_path: Path,
    *,
    left: tuple[int, int],
    right: tuple[int, int],
    extra_args: list[str],
) -> int:
    left_svcf = tmp_path / "left.svcf"
    right_svcf = tmp_path / "right.svcf"

    suffix = "_".join(arg.replace("-", "") for arg in extra_args) or "default"
    output = tmp_path / f"out_{suffix}.svcf"

    _write_del(
        left_svcf,
        record_id="left",
        pos=left[0],
        end=left[1],
        source="callerA",
    )
    _write_del(
        right_svcf,
        record_id="right",
        pos=right[0],
        end=right[1],
        source="callerB",
    )

    _invoke(
        [
            "merge",
            "-i",
            str(left_svcf),
            "-i",
            str(right_svcf),
            "-o",
            str(output),
            "--mode",
            "caller",
            "--caller-names",
            "callerA,callerB",
            "--union",
            *extra_args,
        ]
    )

    return _record_count(output)


def test_issue_187_cli_merge_options_change_output(tmp_path: Path):
    assert (
        _merge_pair(
            tmp_path / "distance_default",
            left=(100, 200),
            right=(300, 400),
            extra_args=[],
        )
        == 2
    )

    assert (
        _merge_pair(
            tmp_path / "distance_override",
            left=(100, 200),
            right=(300, 400),
            extra_args=["--max-distance", "250"],
        )
        == 1
    )

    assert (
        _merge_pair(
            tmp_path / "ratio_default",
            left=(100, 200),
            right=(100, 280),
            extra_args=[],
        )
        == 1
    )

    assert (
        _merge_pair(
            tmp_path / "ratio_override",
            left=(100, 200),
            right=(100, 280),
            extra_args=["--max-length-ratio", "1.5"],
        )
        == 2
    )

    assert (
        _merge_pair(
            tmp_path / "jaccard_default",
            left=(100, 125),
            right=(150, 175),
            extra_args=[],
        )
        == 1
    )

    assert (
        _merge_pair(
            tmp_path / "jaccard_override",
            left=(100, 125),
            right=(150, 175),
            extra_args=["--min-jaccard", "0.05"],
        )
        == 2
    )
