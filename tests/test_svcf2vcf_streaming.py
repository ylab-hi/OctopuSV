from pathlib import Path

import pytest

from octopusv.formatter.svcf_to_vcf_converter import SVCFtoVCFConverter
from octopusv.utils.svcf_parser import SVCFFileEventCreator


FORMAT = "GT:AD:LN:ST:QV:TY:ID:SC:REF:ALT:CO"


def evidence(gt, ad, ln, record_id, source, co="chr1_100-chr1_150"):
    return (
        f"{gt}:{ad}:{ln}:.:60:DEL:{record_id}:{source}:"
        f"N:<DEL>:{co}"
    )


def write_input(path: Path, *, second_row_bad=False):
    header = [
        "##fileformat=VCFv4.2",
        "##source=OctopuSV",
        "##contig=<ID=chr1,length=1000000>",
        "##OctopuSV_mode=multi",
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tA\tB",
    ]
    info1 = (
        "SVTYPE=DEL;END=150;SVLEN=50;CHR2=chr1;SUPPORT=5;"
        "SVMETHOD=OctopuSV;RTID=.;AF=.;STRAND=.;RNAMES=."
    )
    row1 = (
        f"chr1\t100\tr1\tN\t<DEL>\t60\tPASS\t{info1}\t{FORMAT}\t"
        f"{evidence('0/1', '5,7', '50', 'a1', 'c1')}\t"
        f"{evidence('0/0', '12,0', '50', 'b1', 'c1')}"
    )

    info2 = (
        "SVTYPE=DEL;END=350;SVLEN=50;CHR2=chr1;SUPPORT=6;"
        "SVMETHOD=OctopuSV;RTID=.;AF=.;STRAND=.;RNAMES=."
    )
    row2_blocks = [evidence("0/0", "15,0", "50", "a2", "c1")]
    if not second_row_bad:
        row2_blocks.append(evidence("0/1", "6,8", "50", "b2", "c1"))
    row2 = (
        f"chr1\t300\tr2\tN\t<DEL>\t60\tPASS\t{info2}\t{FORMAT}\t"
        + "\t".join(row2_blocks)
    )

    path.write_text("\n".join(header + [row1, row2]) + "\n")


def test_streaming_output_matches_legacy_materialized_api(tmp_path):
    input_path = tmp_path / "input.svcf"
    output_path = tmp_path / "streamed.vcf"
    write_input(input_path)

    creator = SVCFFileEventCreator([str(input_path)])
    creator.parse()
    materialized = SVCFtoVCFConverter(creator.events, input_path).convert()

    streaming = SVCFtoVCFConverter(None, input_path)
    streaming.convert_to_file(output_path)

    assert output_path.read_text() == materialized


def test_convert_to_file_does_not_call_materializing_convert(tmp_path, monkeypatch):
    input_path = tmp_path / "input.svcf"
    output_path = tmp_path / "output.vcf"
    write_input(input_path)

    converter = SVCFtoVCFConverter(None, input_path)

    def forbidden_convert():
        raise AssertionError("convert() materializes the full output")

    monkeypatch.setattr(converter, "convert", forbidden_convert)
    converter.convert_to_file(output_path)

    assert output_path.exists()
    assert output_path.stat().st_size > 0


def test_failure_does_not_leave_partial_new_output(tmp_path):
    input_path = tmp_path / "bad.svcf"
    output_path = tmp_path / "bad.vcf"
    write_input(input_path, second_row_bad=True)

    converter = SVCFtoVCFConverter(None, input_path)
    with pytest.raises(ValueError, match="Sample column count mismatch"):
        converter.convert_to_file(output_path)

    assert not output_path.exists()
    assert not list(tmp_path.glob(".bad.vcf.*.tmp"))


def test_failure_preserves_existing_output(tmp_path):
    input_path = tmp_path / "bad.svcf"
    output_path = tmp_path / "existing.vcf"
    write_input(input_path, second_row_bad=True)
    output_path.write_text("KEEP_ME\n")

    converter = SVCFtoVCFConverter(None, input_path)
    with pytest.raises(ValueError, match="Sample column count mismatch"):
        converter.convert_to_file(output_path)

    assert output_path.read_text() == "KEEP_ME\n"
    assert not list(tmp_path.glob(".existing.vcf.*.tmp"))


def test_success_preserves_existing_output_permissions(tmp_path):
    input_path = tmp_path / "input.svcf"
    output_path = tmp_path / "existing.vcf"
    write_input(input_path)
    output_path.write_text("OLD\n")
    output_path.chmod(0o640)

    SVCFtoVCFConverter(None, input_path).convert_to_file(output_path)

    assert output_path.stat().st_mode & 0o777 == 0o640
    assert output_path.read_text().startswith("##fileformat=VCFv4.2\n")
