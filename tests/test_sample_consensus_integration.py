from pathlib import Path
import subprocess
import sys

from octopusv.utils.svcf_sample_parser import parse_svcf_sample_block


ROOT = Path(__file__).resolve().parents[1]
DATA_DIR = ROOT / "tests" / "data" / "sample_consensus"
INPUT_DIR = DATA_DIR / "input"
BASELINE = DATA_DIR / "baseline_old" / "sample_union_old.svcf"

SAMPLE_FORMAT_V11 = "GT:AD:UC:UV:LN:ST:QV:TY:ID:SC:REF:ALT:CO"


def _run_sample_merge(output_path: Path):
    cmd = [
        sys.executable,
        "-m",
        "octopusv",
        "merge",
        "-i",
        str(INPUT_DIR / "sample_A_callers.svcf"),
        "-i",
        str(INPUT_DIR / "sample_B_callers.svcf"),
        "-o",
        str(output_path),
        "--mode",
        "sample",
        "--sample-names",
        "sample_A,sample_B",
        "--union",
    ]
    result = subprocess.run(
        cmd,
        cwd=ROOT,
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0, result.stderr
    return result


def _records(path: Path):
    rows = {}
    with path.open(encoding="utf-8-sig") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\r\n").split("\t")
            rows[(parts[0], int(parts[1]))] = parts
    return rows


def _parse_info(info_text: str):
    result = {}
    for token in info_text.split(";"):
        if not token:
            continue
        if "=" in token:
            key, value = token.split("=", 1)
            result[key] = value
        else:
            result[token] = True
    return result


def _sample_fields(record, sample_index: int):
    fmt = record[8]
    return parse_svcf_sample_block(fmt, record[9 + sample_index])


def test_sample_mode_consensus_uses_svcf_11_sample_schema(tmp_path):
    output = tmp_path / "population.svcf"
    _run_sample_merge(output)

    records = _records(output)
    assert records
    assert {record[8] for record in records.values()} == {SAMPLE_FORMAT_V11}


def test_sample_mode_consensus_expected_gt_uc_uv_and_no_ad(tmp_path):
    output = tmp_path / "population.svcf"
    _run_sample_merge(output)
    records = _records(output)

    # Same expectation for sample_A and sample_B despite deliberately reversed
    # caller order in several input records.
    expected = {
        ("chr1", 1000): ("0/1", "3", "3"),
        ("chr1", 2000): ("1/1", "3", "3"),
        ("chr1", 3000): ("1/.", "3", "3"),
        ("chr1", 4000): ("0/1", "2", "2"),
        ("chr1", 5000): ("1/.", "2", "2"),
        ("chrY", 6000): ("1", "2", "2"),
        ("chr1", 7000): ("0/1", "1", "1"),
        ("chr1", 8000): ("./.", "1", "2"),
    }

    for key, (gt, uc, uv) in expected.items():
        record = records[key]
        for sample_index in (0, 1):
            sample = _sample_fields(record, sample_index)
            assert sample["GT"] == gt, (key, sample_index, sample)
            assert sample["AD"] == ".,.", (key, sample_index, sample)
            assert sample["UC"] == uc, (key, sample_index, sample)
            assert sample["UV"] == uv, (key, sample_index, sample)
            assert sample["SC"] == "OctopuSV", (key, sample_index, sample)


def test_sample_mode_consensus_uses_sample_record_identity_not_first_caller(tmp_path):
    output = tmp_path / "population.svcf"
    _run_sample_merge(output)
    records = _records(output)

    expected_ids = {
        ("chr1", 1000): ("A.ALL_HET", "B.ALL_HET"),
        ("chr1", 2000): ("A.ALL_HOM", "B.ALL_HOM"),
        ("chr1", 3000): ("A.ZYG_DISCORD", "B.ZYG_DISCORD"),
        ("chr1", 4000): ("A.MISSING_CALLER", "B.MISSING_CALLER"),
        ("chr1", 5000): ("A.DUP_SOURCE", "B.DUP_SOURCE"),
        ("chrY", 6000): ("A.HAPLOID", "B.HAPLOID"),
        ("chr1", 7000): ("A.SINGLE_CALLER", "B.SINGLE_CALLER"),
        ("chr1", 8000): ("A.PRESENCE_TIE", "B.PRESENCE_TIE"),
    }

    for key, ids in expected_ids.items():
        record = records[key]
        sample_a = _sample_fields(record, 0)
        sample_b = _sample_fields(record, 1)
        assert sample_a["ID"] == ids[0]
        assert sample_b["ID"] == ids[1]

        info = _parse_info(record[7])
        assert info["SOURCES"] == "sample_A,sample_B"
        assert info["SOURCE_IDS"] == f"{ids[0]},{ids[1]}"


def test_sample_consensus_changes_only_sample_representation_not_core_event(tmp_path):
    output = tmp_path / "population.svcf"
    _run_sample_merge(output)

    old_records = _records(BASELINE)
    new_records = _records(output)
    assert set(new_records) == set(old_records)

    for key in sorted(old_records):
        old = old_records[key]
        new = new_records[key]

        # Fixed record identity/geometry/quality remain unchanged.
        assert new[:7] == old[:7], key

        old_info = _parse_info(old[7])
        new_info = _parse_info(new[7])
        for source_key in ("SOURCES", "SOURCE_IDS"):
            old_info.pop(source_key, None)
            new_info.pop(source_key, None)
        assert new_info == old_info, key


def test_sample_consensus_retires_first_evidence_order_dependence(tmp_path):
    output = tmp_path / "population.svcf"
    result = _run_sample_merge(output)
    records = _records(output)

    zyg = records[("chr1", 3000)]
    assert _sample_fields(zyg, 0)["GT"] == "1/."
    assert _sample_fields(zyg, 1)["GT"] == "1/."

    # The historical first-evidence collapse warning must disappear once
    # consensus synthesis replaces that behavior.
    assert "first deterministic evidence block" not in result.stderr
    assert "collapsed multiple evidence blocks" not in result.stderr
