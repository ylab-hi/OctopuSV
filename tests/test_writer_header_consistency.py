from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace

from octopusv.merger.multi_sample_writer import MultiSampleWriter
from octopusv.merger.name_mapper import NameMapper
from octopusv.utils.svcf_utils import get_octopus_default_definitions


def _write_input_header(path: Path, *, precise_description: str = "Precise structural variant") -> None:
    path.write_text(
        "\n".join(
            [
                "##fileformat=VCFv4.2",
                "##INFO=<ID=PRECISE,Number=0,Type=Flag,Description=\"%s\">" % precise_description,
                "##FILTER=<ID=LowQual,Description=\"Low quality\">",
                "##ALT=<ID=DEL,Description=\"Deletion\">",
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE",
            ]
        )
        + "\n"
    )


def _event():
    return SimpleNamespace(
        chrom="chr1",
        pos=100,
        sv_id="e1",
        ref="N",
        alt="<DEL>",
        quality="152.7",
        filter="LowQual",
        info={
            "SVTYPE": "DEL",
            "END": "200",
            "SVLEN": "100",
            "CHR2": "chr1",
            "SUPPORT": "5",
            "SVMETHOD": "OctopuSV",
            "RTID": ".",
            "AF": "0.5",
            "STRAND": ".",
            "RNAMES": ".",
            "PRECISE": True,
        },
        ordered_samples=[
            {
                "GT": "0/1",
                "AD": ".,.",
                "UC": "1",
                "UV": "1",
                "LN": "100",
                "ST": ".",
                "QV": "152.7",
                "TY": "DEL",
                "ID": "e1",
                "SC": "OctopuSV",
                "REF": "N",
                "ALT": "<DEL>",
                "CO": "chr1_100-chr1_200",
            }
        ],
    )


def test_octopus_default_qv_header_accepts_fractional_quality():
    definitions = get_octopus_default_definitions()["format_lines"]
    qv = next(line for line in definitions if line.startswith("##FORMAT=<ID=QV,"))
    assert "Type=Float" in qv


def test_multi_writer_uses_canonical_af_and_qv_types(tmp_path):
    input_path = tmp_path / "sample.svcf"
    _write_input_header(input_path)
    mapper = NameMapper([str(input_path)], mode="sample")
    writer = MultiSampleWriter(mapper, input_files=[str(input_path)])

    output = tmp_path / "out.svcf"
    writer.write_results(output, [_event()], {"chr1": "1000000"})
    text = output.read_text()

    assert '##INFO=<ID=AF,Number=1,Type=Float' in text
    assert '##FORMAT=<ID=QV,Number=1,Type=Float' in text
    assert '##INFO=<ID=AF,Number=1,Type=String' not in text
    assert '##FORMAT=<ID=QV,Number=1,Type=Integer' not in text


def test_multi_writer_preserves_nonreserved_info_filter_alt_definitions(tmp_path):
    input_path = tmp_path / "sample.svcf"
    _write_input_header(input_path)
    mapper = NameMapper([str(input_path)], mode="sample")
    writer = MultiSampleWriter(mapper, input_files=[str(input_path)])

    output = tmp_path / "out.svcf"
    writer.write_results(output, [_event()], {"chr1": "1000000"})
    text = output.read_text()

    assert '##INFO=<ID=PRECISE,Number=0,Type=Flag,Description="Precise structural variant">' in text
    assert '##FILTER=<ID=LowQual,Description="Low quality">' in text
    assert '##ALT=<ID=DEL,Description="Deletion">' in text
    assert "\tPRECISE;" in text or ";PRECISE;" in text or ";PRECISE\t" in text


def test_multi_writer_does_not_let_input_override_canonical_info_definition(tmp_path):
    input_path = tmp_path / "sample.svcf"
    input_path.write_text(
        "\n".join(
            [
                "##fileformat=VCFv4.2",
                '##INFO=<ID=AF,Number=1,Type=String,Description="Legacy wrong type">',
                '##INFO=<ID=PRECISE,Number=0,Type=Flag,Description="Precise structural variant">',
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE",
            ]
        )
        + "\n"
    )
    mapper = NameMapper([str(input_path)], mode="sample")
    writer = MultiSampleWriter(mapper, input_files=[str(input_path)])

    output = tmp_path / "out.svcf"
    writer.write_results(output, [_event()], {"chr1": "1000000"})
    text = output.read_text()

    assert '##INFO=<ID=AF,Number=1,Type=Float' in text
    assert 'Legacy wrong type' not in text
    assert text.count("##INFO=<ID=AF,") == 1


def _write_valid_caller_svcf(path: Path, *, record_id: str, pos: int) -> None:
    from octopusv.utils.svcf_schema import CALLER_FORMAT

    end = pos + 100
    lines = [
        "##fileformat=VCFv4.2",
        "##SVCFVersion=1.1",
        "##OctopuSV_mode=caller",
        "##contig=<ID=chr1,length=1000000>",
        '##INFO=<ID=PRECISE,Number=0,Type=Flag,Description="Precise structural variant">',
        '##INFO=<ID=AF,Number=1,Type=String,Description="Legacy wrong type">',
        '##FORMAT=<ID=QV,Number=1,Type=Integer,Description="Legacy wrong type">',
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE",
    ]
    info = ";".join(
        [
            "SVTYPE=DEL",
            f"END={end}",
            "SVLEN=100",
            "CHR2=chr1",
            "SUPPORT=5",
            "SVMETHOD=OctopuSV",
            "RTID=.",
            "AF=0.5",
            "STRAND=.",
            "RNAMES=.",
            "PRECISE",
        ]
    )
    sample = f"0/1:5,5:100:.:152.7:DEL:{record_id}:test:N:<DEL>:chr1_{pos}-chr1_{end}"
    lines.append(
        f"chr1\t{pos}\t{record_id}\tN\t<DEL>\t152.7\tPASS\t{info}\t{CALLER_FORMAT}\t{sample}"
    )
    path.write_text("\n".join(lines) + "\n")


def test_sample_merge_wiring_preserves_custom_header_and_canonical_types(tmp_path):
    from typer.testing import CliRunner

    from octopusv.cli.cli import app

    a = tmp_path / "a.svcf"
    b = tmp_path / "b.svcf"
    _write_valid_caller_svcf(a, record_id="a1", pos=100)
    _write_valid_caller_svcf(b, record_id="b1", pos=100)
    output = tmp_path / "population.svcf"

    result = CliRunner().invoke(
        app,
        [
            "merge",
            "-i", str(a),
            "-i", str(b),
            "-o", str(output),
            "--mode", "sample",
            "--sample-names", "A,B",
            "--union",
        ],
    )
    assert result.exit_code == 0, result.output

    text = output.read_text()
    assert '##INFO=<ID=PRECISE,Number=0,Type=Flag,Description="Precise structural variant">' in text
    assert '##INFO=<ID=AF,Number=1,Type=Float' in text
    assert '##FORMAT=<ID=QV,Number=1,Type=Float' in text
    assert 'Legacy wrong type' not in text
    assert "PRECISE" in next(line for line in text.splitlines() if not line.startswith("#")).split("\t")[7]
