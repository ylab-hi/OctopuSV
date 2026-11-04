"""Guard tests on real VCaP fixtures (Stage 16, 10 callers), written before #192.

Each fixture under tests/data/merge_guard_vcap/<case>/ was cut from the real
caller-mode merge inputs by octopusv_guard.py and verified on Quest to
reproduce its baseline groups exactly (see cut_report.txt).  The code that
produced expected_union.svcf had source sha256
743059a9289b5f69abf16ee288bf36e8a3bc43839fbaf8e95e3d29cd271229bd.

Two kinds of tests:

* test_fixture_output_is_locked: the whole merged output is frozen.  Only
  ##fileDate and the absolute input paths inside
  ##OctopuSV_input_reference audit lines are normalised.
* Named behaviour tests: the specific facts each fixture exists to protect,
  so a failure says what broke.

Tests marked STATUS QUO describe behaviour that issue #192 is expected to
change on purpose; update them in that change and state the reason.
"""

from __future__ import annotations

import re
from pathlib import Path

import pytest
from typer.testing import CliRunner

from octopusv.cli.cli import app


DATA = Path(__file__).parent / "data" / "merge_guard_vcap"
CASES = sorted(p.name for p in DATA.iterdir() if (p / "expected_union.svcf").is_file())
SOURCE_PATH_RE = re.compile(r'Source="[^"]*/([^"/]+)"')

runner = CliRunner()
_cache: dict[str, str] = {}


def _normalise(text: str) -> list[str]:
    return [
        SOURCE_PATH_RE.sub(r'Source="\1"', line)
        for line in text.splitlines()
        if not line.startswith("##fileDate=")
    ]


def _merged_text(case: str, tmp_path_factory) -> str:
    if case not in _cache:
        out = tmp_path_factory.mktemp(case) / "merged.svcf"
        result = runner.invoke(
            app,
            [
                "merge", "--mode", "caller", "--union",
                "--input-list", str(DATA / case / "inputs" / "manifest.tsv"),
                "-o", str(out),
            ],
        )
        assert result.exit_code == 0, result.output
        _cache[case] = out.read_text()
    return _cache[case]


def _records(case: str, tmp_path_factory) -> dict[str, dict]:
    records = {}
    for line in _merged_text(case, tmp_path_factory).splitlines():
        if not line or line.startswith("#"):
            continue
        fields = line.split("\t")
        info = {}
        for item in fields[7].split(";"):
            key, _, value = item.partition("=")
            info[key] = value
        records[fields[2]] = {
            "chrom": fields[0],
            "pos": int(fields[1]),
            "alt": fields[4],
            "info": info,
            "sources": info.get("SOURCES", "").split(","),
            "source_ids": info.get("SOURCE_IDS", "").split(","),
            "n_blocks": len(fields) - 9,
        }
    return records


@pytest.fixture(scope="module")
def merged(tmp_path_factory):
    _cache.clear()
    return lambda case: _records(case, tmp_path_factory)


# ---------------------------------------------------------------------------
# Whole-output lock
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("case", CASES)
def test_fixture_output_is_locked(case, tmp_path_factory):
    expected = (DATA / case / "expected_union.svcf").read_text()
    assert _normalise(_merged_text(case, tmp_path_factory)) == _normalise(expected)


# ---------------------------------------------------------------------------
# Representations of one physical structure that must stay merged
# ---------------------------------------------------------------------------

def test_both_manta_ends_of_one_translocation_merge(merged):
    rec = merged("tra_manta_two_ends")["BND00070006"]
    assert (rec["chrom"], rec["pos"], rec["alt"]) == ("chrM", 16093, "[chr11:49862020[T")
    assert rec["info"]["SVTYPE"] == "TRA"
    assert rec["sources"] == ["delly", "lumpy", "manta", "manta"]
    assert rec["source_ids"] == [
        "BND00070006", "3840_1",
        "MantaBND:0:40429:40438:0:1:0:1", "MantaBND:0:40429:40438:0:1:0:0",
    ]


def test_two_delly_halves_of_one_inversion_merge(merged):
    rec = merged("inv_delly_two_halves")["INV00000108"]
    assert (rec["chrom"], rec["pos"], rec["info"]["END"]) == ("chr1", 6558910, "6559722")
    assert rec["sources"] == ["delly", "delly", "lumpy", "svaba"]
    assert rec["source_ids"] == ["INV00000108", "INV00000109", "8_1", "943033226:1"]


def test_pbsv_split_dup_and_bnd_dup_merge(merged):
    rec = merged("dup_pbsv_split_dup_and_bnd")["DUP00041755"]
    assert rec["info"]["SVTYPE"] == "DUP"
    assert (rec["pos"], rec["info"]["END"]) == (132146176, "132357817")
    assert rec["sources"].count("pbsv") == 2
    assert "pbsv.SPLIT.DUP.10090" in rec["source_ids"]
    assert "pbsv.BND.chr3:132146176-chr3:132357823" in rec["source_ids"]
    assert len(rec["sources"]) == 9


def test_large_dup_with_two_records_from_one_caller_merges(merged):
    rec = merged("large_same_source_1")["DUP00041530"]
    assert (rec["pos"], rec["info"]["END"]) == (117321605, "117904941")
    assert rec["sources"] == [
        "cutesv", "debreak", "delly", "lumpy", "manta", "pbsv", "pbsv", "svim",
    ]


def test_sv3_translocation_merges_eight_callers_with_same_representative(merged):
    rec = merged("tra_sv3_zdhhc7_8_callers")["BND00022166"]
    assert (rec["chrom"], rec["pos"], rec["alt"]) == ("chr17", 75792754, "A[chr16:84992411[")
    assert rec["info"]["CHR2"] == "chr16" and rec["info"]["END"] == "84992411"
    assert rec["sources"] == [
        "cutesv", "delly", "lumpy", "manta", "pbsv", "sniffles2", "svaba", "svim",
    ]


# ---------------------------------------------------------------------------
# RamiSV truth SVs: grouping and representative must not move
# ---------------------------------------------------------------------------

RAMISV_TRUTH = {
    ("ramisv_usp10_zdhhc7_sv1_sv2_sv3", "INV00018420"): (
        "chr16", 84712741, "<INV>", "84992555",
        ["cutesv", "debreak", "delly", "lumpy", "manta", "pbsv", "sniffles2", "svaba", "svim"],
    ),
    ("ramisv_usp10_zdhhc7_sv1_sv2_sv3", "BND00018421"): (
        "chr16", 84986672, "]chr12:122972024]C", "122972024",
        ["cutesv", "delly", "lumpy", "manta", "pbsv", "svaba", "svim"],
    ),
    ("ramisv_usp10_zdhhc7_sv1_sv2_sv3", "BND00022166"): (
        "chr17", 75792754, "A[chr16:84992411[", "84992411",
        ["cutesv", "delly", "lumpy", "manta", "pbsv", "sniffles2", "svaba", "svim"],
    ),
    ("ramisv_kpna1_inv_tra", "INV00032147"): (
        "chr2", 204878529, "<INV>", "204879487",
        ["cutesv", "delly", "lumpy", "pbsv", "sniffles2", "svaba", "svim"],
    ),
    ("ramisv_kpna1_inv_tra", "BND00041599"): (
        "chr3", 122424698, "]chr2:204878835]C", "204878835",
        ["cutesv", "delly", "lumpy", "manta", "pbsv", "svim"],
    ),
    ("ramisv_chd1", "DB16961"): (
        "chr5", 68062107, "<INV>", "98966481", ["debreak", "delly", "svim"],
    ),
}


@pytest.mark.parametrize("case,record_id", sorted(RAMISV_TRUTH))
def test_ramisv_truth_sv_group_and_representative(merged, case, record_id):
    chrom, pos, alt, end, sources = RAMISV_TRUTH[(case, record_id)]
    rec = merged(case)[record_id]
    assert (rec["chrom"], rec["pos"], rec["alt"], rec["info"]["END"]) == (chrom, pos, alt, end)
    assert rec["sources"] == sources


# ---------------------------------------------------------------------------
# Current behaviour kept on purpose (not touched by #192)
# ---------------------------------------------------------------------------

def test_manta_single_inverted_junction_stays_two_bnd_records(merged):
    # The two Manta mates are 3 bp short of mirroring each other exactly, so
    # correct keeps them as BND; merge does not join same-chromosome BND mates.
    # Whether such a pair should become one record is a deferred decision.
    records = merged("bnd_manta_single_inverted_junction")
    for record_id in ("MantaBND:109988:3:4:0:0:0:1", "MantaBND:109988:3:4:0:0:0:0"):
        assert records[record_id]["info"]["SVTYPE"] == "BND"
        assert records[record_id]["sources"] == ["manta"]


# ---------------------------------------------------------------------------
# STATUS QUO: expected to change in #192
# ---------------------------------------------------------------------------

def test_chr5_56bp_and_112bp_deletions_status_quo_before_192(merged):
    # STATUS QUO (#192): six callers each report both a 56 bp and a 112 bp
    # deletion here, yet they form one group (seed Sniffles2 56 bp; the DEL
    # length-ratio limit 2.0 is inclusive).  #192 should split them.
    rec = merged("del_chr5_56_112_a")["cuteSV.DEL.15482"]
    assert (rec["pos"], rec["info"]["SVLEN"]) == (177467978, "55")
    assert rec["n_blocks"] == 13
    for caller in ("cutesv", "debreak", "hitsv", "pbsv", "sniffles2", "svim"):
        assert rec["sources"].count(caller) == 2


def test_chr5_second_56_112_site_status_quo_before_192(merged):
    rec = merged("del_chr5_56_112_b")["pbsv.DEL.17403"]
    assert (rec["pos"], rec["info"]["SVLEN"]) == (179001504, "112")
    assert rec["n_blocks"] == 13


def test_inv_00000658_status_quo_pending_review(merged):
    # Two delly and two lumpy records (two inversion halves each).  Kept as is
    # until the two breakpoint sets are reviewed.
    rec = merged("inv_00000658_to_check")["INV00000658"]
    assert rec["sources"] == ["delly", "delly", "lumpy", "lumpy", "manta"]
    assert rec["info"]["SVLEN"] == "1209"
