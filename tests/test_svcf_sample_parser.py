from __future__ import annotations

import pytest

from octopusv.utils.svcf_sample_parser import parse_svcf_sample_block
from octopusv.utils.svcf_parser import SVCFEvent


FORMAT = "GT:AD:LN:ST:QV:TY:ID:SC:REF:ALT:CO"


def test_parse_normal_ins_block():
    block = (
        "0/1:5,7:73:.:60:INS:cuteSV.INS.22:cuteSV-2.0.3:"
        "A:<INS>:chr1_1202256-chr1_1202329"
    )

    parsed = parse_svcf_sample_block(FORMAT, block)

    assert parsed == {
        "GT": "0/1",
        "AD": "5,7",
        "LN": "73",
        "ST": ".",
        "QV": "60",
        "TY": "INS",
        "ID": "cuteSV.INS.22",
        "SC": "cuteSV-2.0.3",
        "REF": "A",
        "ALT": "<INS>",
        "CO": "chr1_1202256-chr1_1202329",
    }


def test_parse_manta_id_with_colons():
    block = (
        "0/1:8,12:101:.:999:DEL:"
        "MantaDEL:469174:0:1:0:0:0:"
        "Manta_v1.6.0:N:<DEL>:chr1_100-chr1_201"
    )

    parsed = parse_svcf_sample_block(FORMAT, block)

    assert parsed["ID"] == "MantaDEL:469174:0:1:0:0:0"
    assert parsed["SC"] == "Manta_v1.6.0"
    assert parsed["REF"] == "N"
    assert parsed["ALT"] == "<DEL>"
    assert parsed["CO"] == "chr1_100-chr1_201"


@pytest.mark.parametrize(
    "alt",
    [
        "N[chr2:12345[",
        "N]chr2:12345]",
        "[chr2:12345[N",
        "]chr2:12345]N",
    ],
)
def test_parse_bnd_alt_with_colon(alt):
    block = (
        "0/1:.,.:.:+-:60:BND:bnd.1:Sniffles2_2.2:N:"
        f"{alt}:chr1_500-chr2_12345"
    )

    parsed = parse_svcf_sample_block(FORMAT, block)

    assert parsed["ID"] == "bnd.1"
    assert parsed["SC"] == "Sniffles2_2.2"
    assert parsed["REF"] == "N"
    assert parsed["ALT"] == alt
    assert parsed["CO"] == "chr1_500-chr2_12345"


def test_parse_manta_id_and_bnd_alt_together():
    alt = "N]chr2:12345]"
    block = (
        "0/1:.,.:.:+-:60:BND:"
        "MantaBND:469174:0:1:0:0:0:"
        "Manta_v1.6.0:N:"
        f"{alt}:chr1_500-chr2_12345"
    )

    parsed = parse_svcf_sample_block(FORMAT, block)

    assert parsed["ID"] == "MantaBND:469174:0:1:0:0:0"
    assert parsed["SC"] == "Manta_v1.6.0"
    assert parsed["REF"] == "N"
    assert parsed["ALT"] == alt
    assert parsed["CO"] == "chr1_500-chr2_12345"


def test_parse_co_with_underscore_contigs():
    co = "chrY_KI270740v1_random_100-NC_007605_250"
    block = (
        "0/1:3,9:150:.:40:DEL:caller.1:caller:N:<DEL>:"
        f"{co}"
    )

    parsed = parse_svcf_sample_block(FORMAT, block)

    assert parsed["CO"] == co


def test_missing_placeholder_keeps_fixed_schema():
    block = "0/0:.:.:.:.:.:.:.:.:.:."

    parsed = parse_svcf_sample_block(FORMAT, block)

    assert list(parsed) == FORMAT.split(":")
    assert parsed["GT"] == "0/0"
    assert parsed["ID"] == "."
    assert parsed["ALT"] == "."
    assert parsed["CO"] == "."


def test_svcf_event_uses_shared_parser_but_keeps_record_original_id():
    block = (
        "0/1:8,12:101:.:999:DEL:"
        "MantaDEL:469174:0:1:0:0:0:"
        "Manta_v1.6.0:N:<DEL>:chr1_100-chr1_201"
    )

    event = SVCFEvent(
        chrom="chr1",
        pos="100",
        sv_id="representative.1",
        ref="N",
        alt="<DEL>",
        quality="999",
        filter="PASS",
        info=(
            "SVTYPE=DEL;END=201;SVLEN=101;CHR2=chr1;SUPPORT=20;"
            "SVMETHOD=OctopuSV;RTID=.;AF=.;STRAND=.;RNAMES=."
        ),
        format=FORMAT,
        sample=block,
        source_file="input.svcf",
        sample_name="SAMPLE",
    )

    assert event.sample["ID"] == "MantaDEL:469174:0:1:0:0:0"
    assert event.sample["SC"] == "Manta_v1.6.0"
    assert event.sample["REF"] == "N"
    assert event.sample["ALT"] == "<DEL>"
    assert event.sample["CO"] == "chr1_100-chr1_201"
    assert event.sample["original_id"] == "representative.1"


def test_inspector_wrapper_uses_shared_parser():
    from octopusv.inspection.svcf_inspector import parse_format_block

    block = (
        "0/1:8,12:101:.:999:DEL:"
        "MantaDEL:469174:0:1:0:0:0:"
        "Manta_v1.6.0:N:<DEL>:chr1_100-chr1_201"
    )

    parsed = parse_format_block(FORMAT.split(":"), block)

    assert parsed["ID"] == "MantaDEL:469174:0:1:0:0:0"
    assert parsed["SC"] == "Manta_v1.6.0"
    assert parsed["CO"] == "chr1_100-chr1_201"


def test_multi_sample_writer_extracts_full_colon_id_from_raw_block():
    from octopusv.merger.multi_sample_writer import MultiSampleWriter

    block = (
        "0/1:8,12:101:.:999:DEL:"
        "MantaDEL:469174:0:1:0:0:0:"
        "Manta_v1.6.0:N:<DEL>:chr1_100-chr1_201"
    )

    writer = MultiSampleWriter(name_mapper=None)

    assert writer._sample_id_from_data(
        block,
        FORMAT.split(":"),
    ) == "MantaDEL:469174:0:1:0:0:0"


def test_parse_symbolic_alt_with_colon_subtypes():
    block = (
        "0/1:5,7:300:.:60:INS:caller.1:caller_v1:N:"
        "<INS:ME:ALU>:chr1_100-chr1_400"
    )

    parsed = parse_svcf_sample_block(FORMAT, block)

    assert parsed["ID"] == "caller.1"
    assert parsed["SC"] == "caller_v1"
    assert parsed["REF"] == "N"
    assert parsed["ALT"] == "<INS:ME:ALU>"
    assert parsed["CO"] == "chr1_100-chr1_400"


def test_parse_colon_id_and_symbolic_alt_with_colons_together():
    block = (
        "0/1:5,7:300:.:60:INS:"
        "MantaINS:123:0:1:0:0:0:"
        "Manta_v1.6.0:N:<INS:ME:ALU>:chr1_100-chr1_400"
    )

    parsed = parse_svcf_sample_block(FORMAT, block)

    assert parsed["ID"] == "MantaINS:123:0:1:0:0:0"
    assert parsed["SC"] == "Manta_v1.6.0"
    assert parsed["REF"] == "N"
    assert parsed["ALT"] == "<INS:ME:ALU>"
    assert parsed["CO"] == "chr1_100-chr1_400"
