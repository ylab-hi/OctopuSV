from pathlib import Path

import pytest

from octopusv.utils.svcf_validator import SVCFValidator


DATA_DIR = Path(__file__).parent / "data" / "sample_consensus"
INPUT_DIR = DATA_DIR / "input"
BASELINE_DIR = DATA_DIR / "baseline_old"


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
    info = {}
    for item in info_text.split(";"):
        if not item:
            continue
        if "=" in item:
            key, value = item.split("=", 1)
            info[key] = value
    return info


@pytest.mark.parametrize(
    "filename",
    ["sample_A_callers.svcf", "sample_B_callers.svcf"],
)
def test_per_sample_inputs_are_valid_caller_merged_svcf(filename):
    validator = SVCFValidator(str(INPUT_DIR / filename))
    validator.validate()

    assert validator.mode == "caller_merge"
    assert validator.records_total == 8
    assert validator.issues == []


def test_fixture_contains_reversed_caller_order_at_zygosity_disagreement():
    sample_a = _records(INPUT_DIR / "sample_A_callers.svcf")
    sample_b = _records(INPUT_DIR / "sample_B_callers.svcf")

    info_a = _parse_info(sample_a[("chr1", 3000)][7])
    info_b = _parse_info(sample_b[("chr1", 3000)][7])

    assert info_a["SOURCES"] == "cuteSV,svim,pbsv"
    assert info_b["SOURCES"] == "pbsv,svim,cuteSV"

    # Same caller-level facts, different evidence order.
    gts_a = [col.split(":", 1)[0] for col in sample_a[("chr1", 3000)][9:]]
    gts_b = [col.split(":", 1)[0] for col in sample_b[("chr1", 3000)][9:]]
    assert gts_a == ["0/1", "0/1", "1/1"]
    assert gts_b == ["1/1", "0/1", "0/1"]


def test_fixture_contains_same_source_multiple_evidence():
    sample_a = _records(INPUT_DIR / "sample_A_callers.svcf")
    info = _parse_info(sample_a[("chr1", 5000)][7])

    assert info["SOURCES"] == "cuteSV,cuteSV,svim"
    assert len(sample_a[("chr1", 5000)][9:]) == 3


def test_fixture_contains_haploid_alt_call():
    sample_a = _records(INPUT_DIR / "sample_A_callers.svcf")
    sample_b = _records(INPUT_DIR / "sample_B_callers.svcf")

    assert [col.split(":", 1)[0] for col in sample_a[("chrY", 6000)][9:]] == [
        "1",
        "1",
    ]
    assert [col.split(":", 1)[0] for col in sample_b[("chrY", 6000)][9:]] == [
        "1",
        "1",
    ]


def test_historical_baseline_documents_first_evidence_behavior_only():
    """Freeze the old output as evidence, not as the future 1.0 contract."""
    baseline = _records(BASELINE_DIR / "sample_union_old.svcf")

    zyg_discord = baseline[("chr1", 3000)]
    presence_tie = baseline[("chr1", 8000)]

    # Old behavior depended on each per-sample file's first evidence block.
    assert zyg_discord[9].split(":", 1)[0] == "0/1"
    assert zyg_discord[10].split(":", 1)[0] == "1/1"
    assert presence_tie[9].split(":", 1)[0] == "0/1"
    assert presence_tie[10].split(":", 1)[0] == "0/0"
