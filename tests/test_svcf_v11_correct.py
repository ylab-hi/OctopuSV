from __future__ import annotations

from pathlib import Path

from octopusv.sv import SVEvent
from octopusv.utils.svcf_schema import (
    CALLER_FORMAT,
    MODE_CALLER,
    MODE_MULTI,
    parse_identity_from_meta_lines,
)
from octopusv.utils.svcf_utils import generate_sv_header, write_sv_vcf
from octopusv.utils.svcf_validator import SVCFValidator


def _write_raw_vcf_header(path: Path, sample_names: list[str]) -> None:
    path.write_text(
        "\n".join(
            [
                "##fileformat=VCFv4.2",
                "##source=testcaller",
                "##contig=<ID=chr1,length=1000000>",
                "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
                "##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"Allelic depths\">",
                "\t".join(
                    [
                        "#CHROM",
                        "POS",
                        "ID",
                        "REF",
                        "ALT",
                        "QUAL",
                        "FILTER",
                        "INFO",
                        "FORMAT",
                        *sample_names,
                    ]
                ),
            ]
        )
        + "\n"
    )


def _meta_lines(header_lines: list[str]) -> list[str]:
    return [line for line in header_lines if line.startswith("##")]


def _make_event(*, samples: list[str]) -> SVEvent:
    event = SVEvent(
        "chr1",
        100,
        "caller.DEL.1",
        "N",
        "<DEL>",
        "60",
        "PASS",
        "SVTYPE=DEL;END=200;SVLEN=-100;CHR2=chr1;STRAND=+-;RE=10",
        format="GT:AD",
        sample=samples[0],
        samples=samples,
    )
    event.source = "testcaller"
    return event


def test_single_sample_correct_header_declares_v11_caller(tmp_path):
    raw_vcf = tmp_path / "single.vcf"
    _write_raw_vcf_header(raw_vcf, ["S1"])

    header = generate_sv_header(
        ["##contig=<ID=chr1,length=1000000>"],
        str(raw_vcf),
    )
    identity = parse_identity_from_meta_lines(_meta_lines(header))

    assert identity.version == "1.1"
    assert identity.mode == MODE_CALLER
    assert header.count("##SVCFVersion=1.1") == 1
    assert header.count("##OctopuSV_mode=caller") == 1
    assert "##OctopuSV_mode=multi" not in header
    assert header[-1].endswith("\tS1")


def test_single_sample_correct_output_is_valid_v11_caller(tmp_path):
    raw_vcf = tmp_path / "single.vcf"
    output_svcf = tmp_path / "single.svcf"
    _write_raw_vcf_header(raw_vcf, ["S1"])

    write_sv_vcf(
        ["##contig=<ID=chr1,length=1000000>"],
        [_make_event(samples=["0/1:8,4"])],
        output_svcf,
        str(raw_vcf),
    )

    text = output_svcf.read_text()
    assert "##SVCFVersion=1.1\n" in text
    assert "##OctopuSV_mode=caller\n" in text

    record = next(
        line for line in text.splitlines() if line and not line.startswith("#")
    ).split("\t")
    assert record[8] == CALLER_FORMAT
    assert len(record[9:]) == 1

    validator = SVCFValidator(str(output_svcf))
    validator.validate()
    assert [issue for issue in validator.issues if issue.level == "error"] == []


def test_multi_sample_correct_stays_unversioned_legacy_multi(tmp_path):
    raw_vcf = tmp_path / "multi.vcf"
    _write_raw_vcf_header(raw_vcf, ["S1", "S2"])

    header = generate_sv_header(
        ["##contig=<ID=chr1,length=1000000>"],
        str(raw_vcf),
    )
    identity = parse_identity_from_meta_lines(_meta_lines(header))

    # This is intentionally NOT SVCF 1.1 multi: columns are biological
    # samples, but each block still uses caller-evidence FORMAT and no UC/UV
    # sample synthesis has occurred.
    assert identity.version is None
    assert identity.mode == MODE_MULTI
    assert not any(line.startswith("##SVCFVersion=") for line in header)
    assert header.count("##OctopuSV_mode=multi") == 1
    assert header[-1].endswith("\tS1\tS2")


def test_multi_sample_correct_records_remain_caller_format(tmp_path):
    raw_vcf = tmp_path / "multi.vcf"
    output_svcf = tmp_path / "multi.svcf"
    _write_raw_vcf_header(raw_vcf, ["S1", "S2"])

    write_sv_vcf(
        ["##contig=<ID=chr1,length=1000000>"],
        [_make_event(samples=["0/1:8,4", "0/0:12,0"])],
        output_svcf,
        str(raw_vcf),
    )

    lines = output_svcf.read_text().splitlines()
    assert "##OctopuSV_mode=multi" in lines
    assert not any(line.startswith("##SVCFVersion=") for line in lines)

    record = next(line for line in lines if line and not line.startswith("#")).split("\t")
    assert record[8] == CALLER_FORMAT
    assert len(record[9:]) == 2


def test_generate_header_without_input_uses_v11_caller_identity():
    header = generate_sv_header(["##contig=<ID=chr1,length=1000000>"])
    identity = parse_identity_from_meta_lines(_meta_lines(header))

    assert identity.version == "1.1"
    assert identity.mode == MODE_CALLER
    assert header[-1].endswith("\tSample")
