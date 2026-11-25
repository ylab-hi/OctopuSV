"""Representative selection: invalid SUPPORT is unknown, ties never re-parse.

SVCF 1.1 SUPPORT is "." or a non-negative integer.  An invalid value is a
non-blocking validator error, so the file still merges: the value is scored
as unknown (like "."), never reinterpreted, and one warning per source names
it.  The tie breaker reuses the values already scored, so neither an invalid
SUPPORT nor a non-numeric QUAL can crash it.
"""

import logging
from types import SimpleNamespace

import pytest

from octopusv.merger import sv_selector
from octopusv.merger.sv_selector import select_representative_sv


def _sv(sv_id, support, qual="30", caller="sniffles", source=None):
    info = {"SVTYPE": "DEL"}
    if support is not None:
        info["SUPPORT"] = support
    return SimpleNamespace(
        info=info, quality=qual, sample={"SC": caller, "ID": sv_id},
        source_file=source or f"{sv_id}.svcf", sample_name="S", format="GT:ID:SC", sv_id=sv_id,
    )


@pytest.fixture(autouse=True)
def _reset_warnings():
    sv_selector._warned_support_sources.clear()
    yield
    sv_selector._warned_support_sources.clear()


@pytest.mark.parametrize("bad", ["12.0", "NA", "-3", " ", True])
def test_invalid_support_scores_as_unknown_with_one_warning(bad, caplog):
    with caplog.at_level(logging.WARNING):
        rep = select_representative_sv([_sv("y", bad), _sv("x", "1")])
    assert rep.sv_id == "x"                       # invalid counts as 0, not as 12
    assert "Invalid SUPPORT" in caplog.text and "y.svcf" in caplog.text


@pytest.mark.parametrize("missing", [None, ".", ""])
def test_missing_support_scores_as_zero(missing, caplog):
    with caplog.at_level(logging.WARNING):
        rep = select_representative_sv([_sv("x", missing), _sv("y", "1")])
    assert rep.sv_id == "y"
    assert "Invalid SUPPORT" not in caplog.text


def test_invalid_and_missing_support_select_identically():
    a = select_representative_sv([_sv("x", "NA", "40"), _sv("y", ".", "40")])
    b = select_representative_sv([_sv("x", ".", "40"), _sv("y", ".", "40")])
    assert a.sv_id == b.sv_id == "x"


def test_one_warning_per_source(caplog):
    with caplog.at_level(logging.WARNING):
        select_representative_sv([_sv("a", "NA", source="f.svcf"), _sv("b", "1.5", source="f.svcf")])
    assert caplog.text.count("Invalid SUPPORT") == 1


def test_tie_broken_by_support_then_quality():
    # Equal total score is reached when support is capped (>= max_support 30).
    rep = select_representative_sv([_sv("x", "30", "40"), _sv("y", "45", "40")])
    assert rep.sv_id == "y"
    rep = select_representative_sv([_sv("x", "31", "40"), _sv("y", "31", "40")])
    assert rep.sv_id == "x"


def test_non_numeric_qual_does_not_crash_tie_breaker():
    rep = select_representative_sv([_sv("x", "5", "abc"), _sv("y", "5", "abc")])
    assert rep.sv_id == "x"
    rep = select_representative_sv([_sv("x", "5", "50"), _sv("y", "5", "abc")])
    assert rep.sv_id == "x"
