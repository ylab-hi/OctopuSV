import logging
from pathlib import Path
from types import SimpleNamespace

import pytest

from octopusv.formatter.svcf_to_vcf_converter import SVCFtoVCFConverter
from octopusv.merger.multi_sample_writer import MultiSampleWriter
from octopusv.merger.name_mapper import NameMapper


FORMAT = "GT:AD:LN:ST:QV:TY:ID:SC:REF:ALT:CO"
FORMAT_V11 = "GT:AD:UC:UV:LN:ST:QV:TY:ID:SC:REF:ALT:CO"


def evidence(
    gt,
    ad,
    ln,
    record_id,
    source,
    *,
    svtype="DEL",
    ref="N",
    alt="<DEL>",
    co="chr1_100-chr1_150",
):
    return (
        f"{gt}:{ad}:{ln}:.:60:{svtype}:{record_id}:{source}:"
        f"{ref}:{alt}:{co}"
    )


def evidence_v11(
    gt,
    ad,
    uc,
    uv,
    ln,
    record_id,
    *,
    source="OctopuSV",
    svtype="DEL",
    ref="N",
    alt="<DEL>",
    co="chr1_100-chr1_600",
):
    return (
        f"{gt}:{ad}:{uc}:{uv}:{ln}:.:60:{svtype}:{record_id}:{source}:"
        f"{ref}:{alt}:{co}"
    )


def write_svcf(path: Path, sample_names, records, *, multi=True):
    lines = [
        "##fileformat=VCFv4.2",
        "##source=OctopuSV",
        "##contig=<ID=chr1,length=1000000>",
    ]
    if multi:
        lines.append("##OctopuSV_mode=multi")
    lines.append(
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"
        + "\t".join(sample_names)
    )
    lines.extend(records)
    path.write_text("\n".join(lines) + "\n")


def data_lines(path: Path):
    return [
        line.rstrip("\n")
        for line in path.read_text().splitlines()
        if line and not line.startswith("#")
    ]


def _sample_call(
    *,
    gt,
    uc,
    uv,
    ln,
    record_id,
    co="chr1_100-chr1_600",
):
    return {
        "GT": gt,
        "AD": ".,.",
        "UC": str(uc),
        "UV": str(uv),
        "LN": str(ln),
        "ST": ".",
        "QV": "60",
        "TY": "DEL",
        "ID": record_id,
        "SC": "OctopuSV",
        "REF": "N",
        "ALT": "<DEL>",
        "CO": co,
    }


def _write_real_sample_mode_svcf(path: Path, sample_names):
    input_files = [path.parent / f"{name}.caller.svcf" for name in sample_names]
    mapper = NameMapper(
        [str(p) for p in input_files],
        mode="sample",
        custom_names=list(sample_names),
    )
    writer = MultiSampleWriter(mapper)

    calls = [
        _sample_call(
            gt="0/1",
            uc=2,
            uv=3,
            ln=500,
            record_id=f"{sample_names[0]}.event",
        )
    ]
    if len(sample_names) > 1:
        calls.append(
            _sample_call(
                gt="1/.",
                uc=3,
                uv=3,
                ln=500,
                record_id=f"{sample_names[1]}.event",
            )
        )

    event = SimpleNamespace(
        chrom="chr1",
        pos=100,
        sv_id="population.1",
        ref="N",
        alt="<DEL>",
        quality="60",
        filter="PASS",
        info={
            "SVTYPE": "DEL",
            "END": "600",
            "SVLEN": "500",
            "CHR2": "chr1",
            "SUPPORT": "5",
            "SVMETHOD": "OctopuSV",
        },
        ordered_samples=calls,
    )
    writer.write_results(path, [event], {"chr1": 1000000})


def test_multisample_order_and_values_are_preserved(tmp_path):
    input_path = tmp_path / "multi.svcf"
    output_path = tmp_path / "multi.vcf"

    samples = ["sampleA", "sampleB", "sampleC"]
    blocks = [
        evidence("0/1", "5,7", "50", "a1", "callerA"),
        evidence("0/0", "12,0", "50", "b1", "callerB"),
        evidence("./.", ".,.", ".", "c1", "callerC", co="."),
    ]
    info = (
        "SVTYPE=DEL;END=150;SVLEN=50;CHR2=chr1;SUPPORT=7;"
        "SVMETHOD=OctopuSV;RTID=.;AF=.;STRAND=.;RNAMES=."
    )
    record = (
        f"chr1\t100\tmerged.1\tN\t<DEL>\t60\tPASS\t{info}\t{FORMAT}\t"
        + "\t".join(blocks)
    )
    write_svcf(input_path, samples, [record])

    converter = SVCFtoVCFConverter(None, input_path)
    converter.convert_to_file(output_path)

    lines = output_path.read_text().splitlines()
    chrom_header = next(line for line in lines if line.startswith("#CHROM"))
    assert chrom_header.split("\t")[9:] == samples

    fields = next(line for line in lines if not line.startswith("#")).split("\t")
    assert fields[8] == "GT:AD:DP:UC:UV:LN"
    assert fields[9:] == [
        "0/1:5,7:12:.:.:50",
        "0/0:12,0:12:.:.:50",
        "./.:.,.:.:.:.:.",
    ]


def test_real_sample_writer_v11_output_keeps_ln_uc_uv_and_missing_dp(tmp_path):
    input_path = tmp_path / "population.svcf"
    output_path = tmp_path / "population.vcf"
    _write_real_sample_mode_svcf(input_path, ["sampleA", "sampleB"])

    SVCFtoVCFConverter(None, input_path).convert_to_file(output_path)

    lines = output_path.read_text().splitlines()
    assert any(line.startswith("##FORMAT=<ID=UC,") for line in lines)
    assert any(line.startswith("##FORMAT=<ID=UV,") for line in lines)

    fields = next(line for line in lines if not line.startswith("#")).split("\t")
    assert fields[8] == "GT:AD:DP:UC:UV:LN"
    assert fields[9] == "0/1:.,.:.:2:3:500"
    assert fields[10] == "1/.:.,.:.:3:3:500"


@pytest.mark.parametrize("sample_name", ["sampleA"])
def test_single_sample_multi_marker_still_uses_sample_mode(tmp_path, sample_name):
    input_path = tmp_path / "single_population.svcf"
    output_path = tmp_path / "single_population.vcf"
    _write_real_sample_mode_svcf(input_path, [sample_name])

    converter = SVCFtoVCFConverter(None, input_path)
    assert converter.mode == "multi"
    converter.convert_to_file(output_path)

    fields = next(
        line for line in output_path.read_text().splitlines()
        if not line.startswith("#")
    ).split("\t")
    assert fields[8] == "GT:AD:DP:UC:UV:LN"
    assert fields[9] == "0/1:.,.:.:2:3:500"


def test_legacy_multisample_without_marker_warns_and_uses_sample_path(
    tmp_path, caplog
):
    input_path = tmp_path / "legacy_multi.svcf"
    output_path = tmp_path / "legacy_multi.vcf"

    samples = ["sampleA", "sampleB"]
    blocks = [
        evidence("0/1", "5,7", "50", "a1", "callerA"),
        evidence("0/0", "12,0", "50", "b1", "callerB"),
    ]
    info = "SVTYPE=DEL;END=150;SVLEN=50;CHR2=chr1;SUPPORT=7;SVMETHOD=OctopuSV"
    record = (
        f"chr1\t100\tlegacy.1\tN\t<DEL>\t60\tPASS\t{info}\t{FORMAT}\t"
        + "\t".join(blocks)
    )
    write_svcf(input_path, samples, [record], multi=False)

    with caplog.at_level(logging.WARNING):
        converter = SVCFtoVCFConverter(None, input_path)

    assert converter.mode == "legacy_multi"
    assert "treating it as legacy sample-mode input" in caplog.text

    converter.convert_to_file(output_path)
    fields = next(
        line for line in output_path.read_text().splitlines()
        if not line.startswith("#")
    ).split("\t")
    assert fields[8] == "GT:AD:DP:UC:UV:LN"
    assert fields[9:] == [
        "0/1:5,7:12:.:.:50",
        "0/0:12,0:12:.:.:50",
    ]


def test_partial_or_missing_ad_produces_missing_dp(tmp_path):
    input_path = tmp_path / "partial_ad.svcf"
    output_path = tmp_path / "partial_ad.vcf"

    samples = ["sampleA", "sampleB"]
    blocks = [
        evidence_v11("0/1", "5,.", "1", "1", "500", "a1"),
        evidence_v11("0/1", ".,.", "1", "1", "500", "b1"),
    ]
    info = "SVTYPE=DEL;END=600;SVLEN=500;CHR2=chr1;SUPPORT=5;SVMETHOD=OctopuSV"
    record = (
        f"chr1\t100\tpartial.1\tN\t<DEL>\t60\tPASS\t{info}\t{FORMAT_V11}\t"
        + "\t".join(blocks)
    )
    write_svcf(input_path, samples, [record])

    SVCFtoVCFConverter(None, input_path).convert_to_file(output_path)

    fields = next(
        line for line in output_path.read_text().splitlines()
        if not line.startswith("#")
    ).split("\t")
    assert fields[9] == "0/1:5,.:.:1:1:500"
    assert fields[10] == "0/1:.,.:.:1:1:500"


def test_274_samples_remain_274_samples(tmp_path):
    input_path = tmp_path / "cohort274.svcf"
    output_path = tmp_path / "cohort274.vcf"

    samples = [f"S{i:03d}" for i in range(274)]
    blocks = [
        evidence(
            "0/1" if i == 173 else "0/0",
            "8,5" if i == 173 else "13,0",
            "50",
            f"id{i}",
            "callerX",
        )
        for i in range(274)
    ]
    info = (
        "SVTYPE=DEL;END=150;SVLEN=50;CHR2=chr1;SUPPORT=5;"
        "SVMETHOD=OctopuSV;RTID=.;AF=.;STRAND=.;RNAMES=."
    )
    record = (
        f"chr1\t100\tcohort.1\tN\t<DEL>\t60\tPASS\t{info}\t{FORMAT}\t"
        + "\t".join(blocks)
    )
    write_svcf(input_path, samples, [record])

    SVCFtoVCFConverter(None, input_path).convert_to_file(output_path)

    lines = output_path.read_text().splitlines()
    header = next(line for line in lines if line.startswith("#CHROM")).split("\t")
    row = next(line for line in lines if not line.startswith("#")).split("\t")

    assert len(header[9:]) == 274
    assert header[9:] == samples
    assert len(row[9:]) == 274
    assert row[9 + 173].startswith("0/1:8,5:13:.:.:50")


def test_multisample_column_mismatch_fails_loudly(tmp_path):
    input_path = tmp_path / "bad.svcf"

    samples = ["sampleA", "sampleB"]
    only_one_block = evidence("0/1", "5,7", "50", "a1", "callerA")
    info = (
        "SVTYPE=DEL;END=150;SVLEN=50;CHR2=chr1;SUPPORT=7;"
        "SVMETHOD=OctopuSV;RTID=.;AF=.;STRAND=.;RNAMES=."
    )
    record = (
        f"chr1\t100\tbad.1\tN\t<DEL>\t60\tPASS\t{info}\t{FORMAT}\t"
        f"{only_one_block}"
    )
    write_svcf(input_path, samples, [record])

    converter = SVCFtoVCFConverter(None, input_path)
    with pytest.raises(ValueError, match="Sample column count mismatch"):
        converter.convert()


def test_caller_duplicate_evidence_does_not_get_extra_genotype_vote(tmp_path):
    input_path = tmp_path / "caller.svcf"
    output_path = tmp_path / "caller.vcf"

    # Source A contributes two evidence blocks. Evidence-count voting would
    # choose 1/1 (A2 + B1), but unique-source voting uses A1 and B1: the vote is
    # tied and A1 wins by AD support (20 vs 5), so the result must be 0/1.
    blocks = [
        evidence("0/1", "5,20", "50", "A1", "A"),
        evidence("1/1", "1,50", "50", "A2", "A"),
        evidence("1/1", "5,5", "50", "B1", "B"),
    ]
    info = (
        "SVTYPE=DEL;END=150;SVLEN=50;CHR2=chr1;SUPPORT=20;"
        "SVMETHOD=OctopuSV;RTID=.;AF=.;STRAND=.;RNAMES=.;"
        "SOURCES=A,A,B;SOURCE_IDS=A1,A2,B1"
    )
    record = (
        f"chr1\t100\tmerged.dup\tN\t<DEL>\t60\tPASS\t{info}\t{FORMAT}\t"
        + "\t".join(blocks)
    )
    write_svcf(input_path, ["SAMPLE"], [record], multi=False)

    SVCFtoVCFConverter(None, input_path).convert_to_file(output_path)

    fields = next(
        line for line in output_path.read_text().splitlines()
        if not line.startswith("#")
    ).split("\t")
    assert fields[9] == "0/1:5,20:25:50"


def test_colon_containing_later_fields_do_not_shift_gt_ad_ln(tmp_path):
    input_path = tmp_path / "colon.svcf"
    output_path = tmp_path / "colon.vcf"

    block = evidence(
        "0/1",
        "7,9",
        "73",
        "MantaDEL:469174:0:1:0:0:0",
        "Manta",
        svtype="BND",
        alt="N]chr2:12345]",
        co="chr1_100-chr2_12345",
    )
    info = (
        "SVTYPE=BND;END=12345;SVLEN=.;CHR2=chr2;SUPPORT=9;"
        "SVMETHOD=OctopuSV;RTID=.;AF=.;STRAND=.;RNAMES=."
    )
    record = (
        f"chr1\t100\tbnd.1\tN\tN]chr2:12345]\t60\tPASS\t{info}\t{FORMAT}\t{block}"
    )
    write_svcf(input_path, ["SAMPLE"], [record], multi=False)

    SVCFtoVCFConverter(None, input_path).convert_to_file(output_path)

    fields = next(
        line for line in output_path.read_text().splitlines()
        if not line.startswith("#")
    ).split("\t")
    assert fields[9] == "0/1:7,9:16:73"


def test_caller_collapse_ad_ln_come_from_unique_source_winning_block(tmp_path):
    input_path = tmp_path / "caller_duplicate_same_gt.svcf"
    output_path = tmp_path / "caller_duplicate_same_gt.vcf"

    # Unique-source voting uses A1 and B1 only. Their GTs tie, so B1 wins by
    # AD support. A2 is a duplicate-source block that shares B1's winning GT
    # and appears earlier in raw evidence order; it must NOT provide AD/LN.
    blocks = [
        evidence("0/1", "20,4", "50", "A1", "A"),
        evidence("1/1", "0,19", "99", "A2", "A"),
        evidence("1/1", "4,25", "120", "B1", "B"),
    ]
    info = (
        "SVTYPE=DEL;END=220;SVLEN=120;CHR2=chr1;SUPPORT=25;"
        "SVMETHOD=OctopuSV;RTID=.;AF=.;STRAND=.;RNAMES=.;"
        "SOURCES=A,A,B;SOURCE_IDS=A1,A2,B1"
    )
    record = (
        f"chr1\t100\tmerged.samegt\tN\t<DEL>\t60\tPASS\t{info}\t{FORMAT}\t"
        + "\t".join(blocks)
    )
    write_svcf(input_path, ["SAMPLE"], [record], multi=False)

    SVCFtoVCFConverter(None, input_path).convert_to_file(output_path)

    fields = next(
        line for line in output_path.read_text().splitlines()
        if not line.startswith("#")
    ).split("\t")
    assert fields[9] == "1/1:4,25:29:120"
