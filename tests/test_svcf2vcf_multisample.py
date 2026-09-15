from pathlib import Path

import pytest

from octopusv.formatter.svcf_to_vcf_converter import SVCFtoVCFConverter


FORMAT = "GT:AD:LN:ST:QV:TY:ID:SC:REF:ALT:CO"


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
    assert fields[8] == "GT:AD:DP:LN"
    assert fields[9:] == [
        "0/1:5,7:12:50",
        "0/0:12,0:12:50",
        "./.:.,.:0:.",
    ]


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
    assert row[9 + 173].startswith("0/1:8,5:13:50")


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

    # Source A contributes two evidence blocks.  Evidence-count voting would
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

    # Unique-source voting uses A1 and B1 only.  Their GTs tie, so B1 wins by
    # AD support.  A2 is a duplicate-source block that shares B1's winning GT
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
