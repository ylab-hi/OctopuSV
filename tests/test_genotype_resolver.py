from types import SimpleNamespace

from octopusv.formatter.svcf_to_vcf_converter import SVCFtoVCFConverter
from octopusv.utils.genotype_resolver import (
    resolve_multi_caller_genotype,
    unique_source_segments,
)


FORMAT = "GT:AD:LN:ST:QV:TY:ID:SC:REF:ALT:CO"


def _segment(gt, ad, ln, record_id, caller):
    return (
        f"{gt}:{ad}:{ln}:.:60:INS:{record_id}:{caller}:"
        f"N:<INS>:chr1_100-chr1_200"
    )


def test_unique_source_segments_keeps_first_block_per_source():
    segments = ["A1", "B1", "B2", "C1"]
    info = "SOURCES=callerA,callerB,callerB,callerC"

    assert unique_source_segments(segments, info) == [
        (0, "A1"),
        (1, "B1"),
        (3, "C1"),
    ]


def test_duplicate_source_does_not_receive_extra_genotype_vote():
    segments = [
        _segment("1/1", "0,44", 70, "sniffles.1", "Sniffles2"),
        _segment("0/1", "36,9", 140, "svim.1", "SVIM"),
        _segment("0/1", "20,10", 70, "svim.2", "SVIM"),
    ]
    info = "SOURCES=sniffles,svim,svim"

    # Unique callers vote 1/1 vs 0/1. AD breaks the tie in favor of 1/1.
    assert resolve_multi_caller_genotype(FORMAT, segments, info) == "1/1"


def test_unique_source_inputs_keep_existing_majority_rule():
    segments = [
        _segment("0/1", "10,8", 100, "a", "A"),
        _segment("1/1", "0,20", 100, "b", "B"),
        _segment("1/1", "0,15", 100, "c", "C"),
    ]
    info = "SOURCES=A,B,C"

    assert resolve_multi_caller_genotype(FORMAT, segments, info) == "1/1"


def test_missing_or_mismatched_sources_does_not_guess():
    segments = [
        _segment("0/1", "10,8", 100, "a", "A"),
        _segment("1/1", "0,20", 100, "b", "B"),
    ]

    assert unique_source_segments(segments, "SOURCES=A") == [
        (0, segments[0]),
        (1, segments[1]),
    ]


def test_svcf2vcf_collapse_uses_shared_order_independent_consensus():
    converter = object.__new__(SVCFtoVCFConverter)

    segments = [
        _segment("1/1", "0,44", 70, "sniffles.1", "Sniffles2"),
        _segment("0/1", "36,9", 140, "svim.1", "SVIM"),
        _segment("0/1", "20,10", 70, "svim.2", "SVIM"),
    ]

    event = SimpleNamespace(
        format=FORMAT,
        info={
            "SOURCES": "sniffles,svim,svim",
            "SVLEN": "100",
        },
        sv_type="INS",
        pos=100,
        end_pos=100,
    )

    # Sniffles contributes HOM-alt and SVIM contributes one reduced HET state.
    # Presence is established, but zygosity disagrees, so the shared consensus
    # returns 1/. rather than letting AD or input order choose 1/1 or 0/1.
    assert converter._collapse_caller_blocks(event, segments) == "1/.:.,.:.:2:2:100"
