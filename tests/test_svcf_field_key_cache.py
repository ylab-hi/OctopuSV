"""Protect fixed SVCF field-name reuse without changing parsed data."""

from octopusv.utils.svcf_parser import SVCFEvent, parse_info_field
from octopusv.utils.svcf_sample_parser import parse_svcf_sample_block


CALLER_FORMAT = "GT:AD:LN:ST:QV:TY:ID:SC:REF:ALT:CO"


def test_info_fields_keep_exact_values_duplicate_and_unknown_keys():
    raw = (
        "SVTYPE=DEL;END=130;SVLEN=30;CHR2=chr1;SUPPORT=5;"
        "CUSTOM_THING=a=b;FLAG;SVLEN=31;STRAND=."
    )
    parsed = parse_info_field(raw)
    assert parsed == {
        "SVTYPE": "DEL", "END": "130", "SVLEN": "31", "CHR2": "chr1",
        "SUPPORT": "5", "CUSTOM_THING": "a=b", "FLAG": True,
        "STRAND": ".",
    }
    second = parse_info_field(raw)
    first_key = next(key for key in parsed if key == "SVTYPE")
    second_key = next(key for key in second if key == "SVTYPE")
    assert first_key is second_key


def test_format_parsing_preserves_colons_in_source_ids_and_bnd_alt():
    block = (
        "0/1:4,9:12:+-:80:BND:manta:caller:10:GRIDSS2:"
        "A:N]chr2:2345]:chr1_100-chr2_2345"
    )
    expected = {
        "GT": "0/1", "AD": "4,9", "LN": "12", "ST": "+-",
        "QV": "80", "TY": "BND", "ID": "manta:caller:10",
        "SC": "GRIDSS2", "REF": "A", "ALT": "N]chr2:2345]",
        "CO": "chr1_100-chr2_2345",
    }
    first = parse_svcf_sample_block(CALLER_FORMAT, block)
    second = parse_svcf_sample_block(CALLER_FORMAT, block)
    assert first == expected
    assert second == expected
    assert next(key for key in first if key == "GT") is next(key for key in second if key == "GT")


def test_caller_specific_format_does_not_leak_between_schemas():
    block = "1:foo:bar"
    assert parse_svcf_sample_block("GT:AA:BB", block) == {
        "GT": "1", "AA": "foo", "BB": "bar",
    }
    assert parse_svcf_sample_block("GT:BB:AA", block) == {
        "GT": "1", "BB": "foo", "AA": "bar",
    }
    assert parse_svcf_sample_block(["GT", "AA", "BB"], block) == {
        "GT": "1", "AA": "foo", "BB": "bar",
    }


def test_event_parses_coordinates_and_sample_as_before():
    event = SVCFEvent(
        "chr1", "100", "id:1", "A", "<DEL>", "55", "PASS",
        "SVTYPE=DEL;END=130;SVLEN=30;CHR2=chr1;SUPPORT=3",
        CALLER_FORMAT,
        "0/1:4,8:30:.:55:DEL:id:1:cuteSV:A:<DEL>:chr1_100-chr1_130",
        "/data/cutesv.svcf", "SAMPLE",
    )
    assert (event.start_chrom, event.start_pos, event.end_chrom, event.end_pos) == ("chr1", 100, "chr1", 130)
    assert event.sample["ID"] == "id:1"
    assert event.info["SUPPORT"] == "3"
