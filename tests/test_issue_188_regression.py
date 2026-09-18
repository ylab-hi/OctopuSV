from pathlib import Path

from octopusv.formatter.svcf_to_vcf_converter import SVCFtoVCFConverter


CALLER_FORMAT = "GT:AD:LN:ST:QV:TY:ID:SC:REF:ALT:CO"


def _evidence(*, record_id: str, caller: str, co: str) -> str:
    return ":".join(
        [
            "0/1",
            "5,7",
            "26",
            ".",
            "60",
            "DEL",
            record_id,
            caller,
            "N",
            "<DEL>",
            co,
        ]
    )


def _parse_info(info_text: str) -> dict[str, str]:
    parsed = {}
    for item in info_text.split(";"):
        if "=" in item:
            key, value = item.split("=", 1)
            parsed[key] = value
    return parsed


def test_issue_188_svcf2vcf_uses_record_end_not_evidence_co(tmp_path: Path):
    """Merged record END must not be replaced by another caller's CO endpoint."""
    input_svcf = tmp_path / "issue188.svcf"
    output_vcf = tmp_path / "issue188.vcf"

    representative_pos = 125_179_900
    representative_end = 125_179_926

    first_evidence = _evidence(
        record_id="callerA.del",
        caller="callerA",
        co="chr1_125179790-chr1_125179816",
    )
    second_evidence = _evidence(
        record_id="callerB.del",
        caller="callerB",
        co="chr1_125179900-chr1_125179926",
    )

    info = ";".join(
        [
            "SVTYPE=DEL",
            f"END={representative_end}",
            "SVLEN=26",
            "CHR2=chr1",
            "SUPPORT=12",
            "SVMETHOD=OctopuSV",
            "RTID=.",
            "AF=.",
            "STRAND=.",
            "RNAMES=.",
            "SOURCES=callerA,callerB",
            "SOURCE_IDS=callerA.del,callerB.del",
        ]
    )

    input_svcf.write_text(
        "\n".join(
            [
                "##fileformat=VCFv4.2",
                "##SVCFVersion=1.1",
                "##OctopuSV_mode=caller",
                "##contig=<ID=chr1,length=248956422>",
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE",
                (
                    f"chr1\t{representative_pos}\tmerged.del\tN\t<DEL>\t60\tPASS\t"
                    f"{info}\t{CALLER_FORMAT}\t{first_evidence}\t{second_evidence}"
                ),
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    SVCFtoVCFConverter(None, input_svcf).convert_to_file(output_vcf)

    rows = [
        line.split("\t")
        for line in output_vcf.read_text(encoding="utf-8").splitlines()
        if line and not line.startswith("#")
    ]
    assert len(rows) == 1

    out_info = _parse_info(rows[0][7])

    assert rows[0][1] == str(representative_pos)
    assert out_info["END"] == str(representative_end)
    assert out_info["END"] != "125179816"
