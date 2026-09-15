"""SVCF 1.1 sample-mode FORMAT compatibility tests.

UC/UV are intentionally inserted before LN so the fixed evidence tail remains
ID:SC:REF:ALT:CO.  The existing shared parser must therefore retain its robust
handling of colon-containing IDs and ALT values without modification.
"""

from __future__ import annotations

from octopusv.utils.svcf_sample_parser import parse_svcf_sample_block


SAMPLE_FORMAT_V11 = "GT:AD:UC:UV:LN:ST:QV:TY:ID:SC:REF:ALT:CO"


def test_v11_sample_format_preserves_manta_colon_id():
    block = (
        "1/.:.,.:2:3:101:.:999:DEL:"
        "MantaDEL:469174:0:1:0:0:0:"
        "OctopuSV:N:<DEL>:chr1_100-chr1_201"
    )

    parsed = parse_svcf_sample_block(SAMPLE_FORMAT_V11, block)

    assert list(parsed) == SAMPLE_FORMAT_V11.split(":")
    assert parsed["GT"] == "1/."
    assert parsed["AD"] == ".,."
    assert parsed["UC"] == "2"
    assert parsed["UV"] == "3"
    assert parsed["LN"] == "101"
    assert parsed["ID"] == "MantaDEL:469174:0:1:0:0:0"
    assert parsed["SC"] == "OctopuSV"
    assert parsed["REF"] == "N"
    assert parsed["ALT"] == "<DEL>"
    assert parsed["CO"] == "chr1_100-chr1_201"


def test_v11_sample_format_preserves_bnd_alt_with_colon():
    alt = "N]chr2:12345]"
    block = (
        "1/.:.,.:2:2:.:+-:60:BND:"
        "MantaBND:123:0:1:0:0:0:"
        "OctopuSV:N:"
        f"{alt}:chr1_500-chr2_12345"
    )

    parsed = parse_svcf_sample_block(SAMPLE_FORMAT_V11, block)

    assert parsed["GT"] == "1/."
    assert parsed["UC"] == "2"
    assert parsed["UV"] == "2"
    assert parsed["ID"] == "MantaBND:123:0:1:0:0:0"
    assert parsed["SC"] == "OctopuSV"
    assert parsed["REF"] == "N"
    assert parsed["ALT"] == alt
    assert parsed["CO"] == "chr1_500-chr2_12345"
