from __future__ import annotations

from collections import Counter
import os
import subprocess
import sys
from pathlib import Path

import pytest
from typer.testing import CliRunner

pytest.importorskip("natsort")

from octopusv.cli.cli import app
from octopusv.utils.svcf_sample_parser import parse_svcf_sample_block


RUNNER = CliRunner()
CALLER_FORMAT = "GT:AD:LN:ST:QV:TY:ID:SC:REF:ALT:CO"
REPO_ROOT = Path(__file__).resolve().parents[1]


def _invoke(args: list[str]):
    result = RUNNER.invoke(app, args)
    assert result.exit_code == 0, (
        f"Command failed: octopusv {' '.join(args)}\n"
        f"output:\n{result.output}\n"
        f"exception: {result.exception!r}"
    )
    return result


def _parse_info(text: str) -> dict[str, str | bool]:
    info: dict[str, str | bool] = {}
    for item in text.split(";"):
        if "=" in item:
            key, value = item.split("=", 1)
            info[key] = value
        elif item:
            info[item] = True
    return info


def _data_rows(path: Path) -> list[list[str]]:
    return [
        line.split("\t")
        for line in path.read_text(encoding="utf-8").splitlines()
        if line and not line.startswith("#")
    ]


def _bnd(chrom: str, pos: int, record_id: str, alt: str) -> str:
    return (
        f"{chrom}\t{pos}\t{record_id}\tN\t{alt}\t60\tPASS\t"
        "SVTYPE=BND\tGT\t0/1\n"
    )


def _raw_sv(
    chrom: str,
    pos: int,
    record_id: str,
    svtype: str,
    *,
    end: int,
    svlen: int,
) -> str:
    return (
        f"{chrom}\t{pos}\t{record_id}\tN\t<{svtype}>\t60\tPASS\t"
        f"SVTYPE={svtype};END={end};SVLEN={svlen}\tGT\t0/1\n"
    )


def _write_raw_vcf(path: Path, records: list[str]) -> None:
    path.write_text(
        "##fileformat=VCFv4.2\n"
        "##source=MaintenanceInvariantCaller\n"
        "##contig=<ID=chr1,length=1000000>\n"
        "##contig=<ID=chr2,length=1000000>\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n"
        + "".join(records),
        encoding="utf-8",
    )


def _correct_semantic_signature(path: Path):
    """Compare corrected biological records, not their original input order."""
    signatures = []
    for row in _data_rows(path):
        info = _parse_info(row[7])
        signatures.append(
            (
                row[0],
                row[1],
                row[2],
                row[3],
                row[4],
                str(info["SVTYPE"]),
                str(info["END"]),
                str(info["SVLEN"]),
                str(info["CHR2"]),
                str(info["RTID"]),
            )
        )
    return Counter(signatures)


def test_correct_raw_record_order_does_not_change_unambiguous_event_semantics(tmp_path):
    """Changing raw-record order must not change clearly defined corrected events.

    The fixture deliberately avoids the separately documented ambiguous/pending
    BND-pair characterization.  It mixes native SVs, safe same-chromosome BND
    conversions, a single interchromosomal BND->TRA, and a retained BND.
    """
    records = [
        _raw_sv("chr1", 1000, "native.del", "DEL", end=1100, svlen=-100),
        _raw_sv("chr1", 2000, "native.ins", "INS", end=2000, svlen=40),
        _bnd("chr1", 3000, "pair.del.left", "N[chr1:3100["),
        _bnd("chr1", 3100, "pair.del.right", "]chr1:3000]N"),
        _bnd("chr1", 5000, "pair.inv.left", "N]chr1:5200]"),
        _bnd("chr1", 5200, "pair.inv.right", "N]chr1:5000]"),
        _bnd("chr1", 7000, "single.tra", "N[chr2:8000["),
        _bnd("chr1", 9000, "kept.bnd", "N[chr1:9300["),
    ]
    orders = [
        records,
        list(reversed(records)),
        [records[index] for index in (6, 2, 0, 7, 4, 1, 3, 5)],
    ]

    signatures = []
    for index, ordered_records in enumerate(orders):
        raw = tmp_path / f"raw_{index}.vcf"
        corrected = tmp_path / f"corrected_{index}.svcf"
        _write_raw_vcf(raw, ordered_records)
        _invoke(["correct", "-i", str(raw), "-o", str(corrected)])
        _invoke(["validate-svcf", str(corrected)])
        signatures.append(_correct_semantic_signature(corrected))

    assert signatures[0] == signatures[1] == signatures[2]


def _write_caller_svcf(
    path: Path,
    *,
    caller: str,
    records: list[tuple[str, int, int]],
) -> None:
    lines = [
        "##fileformat=VCFv4.2\n",
        "##SVCFVersion=1.1\n",
        "##OctopuSV_mode=caller\n",
        "##contig=<ID=chr1,length=1000000>\n",
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n",
    ]
    for record_id, pos, end in records:
        length = end - pos
        info = (
            f"SVTYPE=DEL;END={end};SVLEN={length};CHR2=chr1;SUPPORT=5;"
            "SVMETHOD=OctopuSV;RTID=.;AF=.;STRAND=.;RNAMES=."
        )
        block = (
            f"0/1:5,5:{length}:.:60:DEL:{record_id}:{caller}:N:<DEL>:"
            f"chr1_{pos}-chr1_{end}"
        )
        lines.append(
            f"chr1\t{pos}\t{record_id}\tN\t<DEL>\t60\tPASS\t{info}\t"
            f"{CALLER_FORMAT}\t{block}\n"
        )
    path.write_text("".join(lines), encoding="utf-8")


def _build_three_inputs(tmp_path: Path) -> tuple[Path, Path, Path]:
    source_a = tmp_path / "a.svcf"
    source_b = tmp_path / "b.svcf"
    source_c = tmp_path / "c.svcf"
    _write_caller_svcf(
        source_a,
        caller="callerA",
        records=[
            ("A.shared", 100, 200),
            ("A.ab", 1000, 1100),
            ("A.only", 2000, 2100),
        ],
    )
    _write_caller_svcf(
        source_b,
        caller="callerB",
        records=[
            ("B.shared", 105, 205),
            ("B.ab", 1005, 1105),
            ("B.only", 3000, 3100),
        ],
    )
    _write_caller_svcf(
        source_c,
        caller="callerC",
        records=[
            ("C.shared", 110, 210),
            ("C.only", 4000, 4100),
        ],
    )
    return source_a, source_b, source_c


def _source_ids(row: list[str]) -> tuple[str, ...]:
    info = _parse_info(row[7])
    return tuple(str(info["SOURCE_IDS"]).split(","))


def _event_geometry(row: list[str]):
    info = _parse_info(row[7])
    return (
        row[0],
        row[1],
        row[3],
        row[4],
        str(info["SVTYPE"]),
        str(info["END"]),
        str(info["SVLEN"]),
        str(info["CHR2"]),
    )


def _membership_geometry_map(path: Path):
    """Map exact input-event membership to merged-event geometry."""
    return {
        frozenset(_source_ids(row)): _event_geometry(row)
        for row in _data_rows(path)
    }


def _caller_group_content_map(path: Path):
    """Capture group content/binding while ignoring unrelated header details."""
    result = {}
    for row in _data_rows(path):
        info = _parse_info(row[7])
        source_ids = tuple(str(info["SOURCE_IDS"]).split(","))
        sources = tuple(str(info["SOURCES"]).split(","))
        evidence = tuple(
            (
                str(block["ID"]),
                str(block["SC"]),
                str(block["GT"]),
                str(block["TY"]),
            )
            for block in (
                parse_svcf_sample_block(row[8], raw_block)
                for raw_block in row[9:]
            )
        )
        result[frozenset(source_ids)] = (
            _event_geometry(row),
            sources,
            source_ids,
            evidence,
        )
    return result


def test_caller_and_sample_modes_share_the_same_event_identity_layer(tmp_path):
    """Mode may change representation, but not which input events form a group."""
    source_a, source_b, source_c = _build_three_inputs(tmp_path)
    inputs = [source_a, source_b, source_c]

    caller_out = tmp_path / "caller_union.svcf"
    _invoke(
        [
            "merge",
            *sum((["-i", str(path)] for path in inputs), []),
            "-o", str(caller_out),
            "--mode", "caller",
            "--caller-names", "A,B,C",
            "--union",
        ]
    )

    sample_out = tmp_path / "sample_union.svcf"
    _invoke(
        [
            "merge",
            *sum((["-i", str(path)] for path in inputs), []),
            "-o", str(sample_out),
            "--mode", "sample",
            "--sample-names", "A,B,C",
            "--union",
        ]
    )

    _invoke(["validate-svcf", str(caller_out)])
    _invoke(["validate-svcf", str(sample_out)])

    assert _membership_geometry_map(caller_out) == _membership_geometry_map(sample_out)


def test_merge_selection_changes_the_event_set_not_surviving_group_contents(tmp_path):
    """Selection acts after grouping: surviving groups must remain unchanged."""
    source_a, source_b, source_c = _build_three_inputs(tmp_path)
    inputs = [source_a, source_b, source_c]
    input_args = sum((["-i", str(path)] for path in inputs), [])

    union_out = tmp_path / "union.svcf"
    _invoke(
        [
            "merge", *input_args,
            "-o", str(union_out),
            "--mode", "caller",
            "--caller-names", "A,B,C",
            "--union",
        ]
    )
    union_groups = _caller_group_content_map(union_out)

    cases = [
        ("min2", ["--min-support", "2"], {
            frozenset({"A.shared", "B.shared", "C.shared"}),
            frozenset({"A.ab", "B.ab"}),
        }),
        ("intersection", ["--intersect"], {
            frozenset({"A.shared", "B.shared", "C.shared"}),
        }),
        ("exact1", ["--exact-support", "1"], {
            frozenset({"A.only"}),
            frozenset({"B.only"}),
            frozenset({"C.only"}),
        }),
    ]

    for name, selection_args, expected_memberships in cases:
        selected_out = tmp_path / f"{name}.svcf"
        _invoke(
            [
                "merge", *input_args,
                "-o", str(selected_out),
                "--mode", "caller",
                "--caller-names", "A,B,C",
                *selection_args,
            ]
        )
        selected_groups = _caller_group_content_map(selected_out)
        assert set(selected_groups) == expected_memberships
        for membership, content in selected_groups.items():
            assert content == union_groups[membership]


def _label_independent_group_map(path: Path):
    """Ignore display SOURCES labels while preserving scientific/binding content."""
    result = {}
    for row in _data_rows(path):
        source_ids = _source_ids(row)
        evidence = tuple(
            (
                str(block["ID"]),
                str(block["SC"]),
                str(block["GT"]),
            )
            for block in (
                parse_svcf_sample_block(row[8], raw_block)
                for raw_block in row[9:]
            )
        )
        result[frozenset(source_ids)] = (_event_geometry(row), source_ids, evidence)
    return result


def test_custom_caller_labels_do_not_participate_in_scientific_grouping(tmp_path):
    source_a, source_b, source_c = _build_three_inputs(tmp_path)
    inputs = [source_a, source_b, source_c]
    input_args = sum((["-i", str(path)] for path in inputs), [])

    first = tmp_path / "labels_abc.svcf"
    _invoke(
        [
            "merge", *input_args,
            "-o", str(first),
            "--mode", "caller",
            "--caller-names", "A,B,C",
            "--union",
        ]
    )

    second = tmp_path / "labels_xyz.svcf"
    _invoke(
        [
            "merge", *input_args,
            "-o", str(second),
            "--mode", "caller",
            "--caller-names", "alpha,beta,gamma",
            "--union",
        ]
    )

    assert _label_independent_group_map(first) == _label_independent_group_map(second)

    first_sources = {
        frozenset(_source_ids(row)): tuple(str(_parse_info(row[7])["SOURCES"]).split(","))
        for row in _data_rows(first)
    }
    second_sources = {
        frozenset(_source_ids(row)): tuple(str(_parse_info(row[7])["SOURCES"]).split(","))
        for row in _data_rows(second)
    }
    assert first_sources != second_sources


def _subprocess_env(hash_seed: str) -> dict[str, str]:
    env = os.environ.copy()
    src = REPO_ROOT / "src"
    package_root = src if src.exists() else REPO_ROOT
    current = env.get("PYTHONPATH", "")
    env["PYTHONPATH"] = str(package_root) + (os.pathsep + current if current else "")
    env["PYTHONHASHSEED"] = hash_seed
    return env


def _run_octopusv_subprocess(args: list[str], *, hash_seed: str) -> None:
    result = subprocess.run(
        [sys.executable, "-m", "octopusv", *args],
        cwd=REPO_ROOT,
        env=_subprocess_env(hash_seed),
        capture_output=True,
        text=True,
        check=False,
    )
    assert result.returncode == 0, (
        f"PYTHONHASHSEED={hash_seed} failed: octopusv {' '.join(args)}\n"
        f"stdout:\n{result.stdout}\n"
        f"stderr:\n{result.stderr}"
    )


def _stable_output_bytes(path: Path) -> bytes:
    """Ignore intentionally volatile run-time metadata, but nothing scientific."""
    lines = path.read_bytes().splitlines(keepends=True)
    return b"".join(line for line in lines if not line.startswith(b"##fileDate="))


def test_merge_and_vcf_export_are_cross_process_deterministic(tmp_path):
    """Unordered Python containers must not make semantic output process-dependent."""
    source_a, source_b, source_c = _build_three_inputs(tmp_path)
    inputs = [source_a, source_b, source_c]

    outputs = []
    for seed in ("1", "777"):
        merged = tmp_path / f"merged_seed_{seed}.svcf"
        final_vcf = tmp_path / f"final_seed_{seed}.vcf"
        merge_args = ["merge"]
        for path in inputs:
            merge_args.extend(["-i", str(path)])
        merge_args.extend(
            [
                "-o", str(merged),
                "--mode", "caller",
                "--caller-names", "A,B,C",
                "--union",
            ]
        )
        _run_octopusv_subprocess(merge_args, hash_seed=seed)
        _run_octopusv_subprocess(
            ["svcf2vcf", "-i", str(merged), "-o", str(final_vcf)],
            hash_seed=seed,
        )
        outputs.append(
            (_stable_output_bytes(merged), _stable_output_bytes(final_vcf))
        )

    assert outputs[0] == outputs[1]
