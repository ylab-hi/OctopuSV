from __future__ import annotations

from collections import Counter
from pathlib import Path
from types import SimpleNamespace

import pytest
from typer.testing import CliRunner

pytest.importorskip("natsort")

from octopusv.cli.cli import app
from octopusv.merger.multi_sample_writer import MultiSampleWriter
from octopusv.merger.name_mapper import NameMapper
from octopusv.merger.sv_merger import SVMerger
from octopusv.sv import SVEvent
from octopusv.utils.svcf_utils import write_sv_vcf


RUNNER = CliRunner()
CALLER_FORMAT = "GT:ID:SC"


def _data_rows(path: Path) -> list[str]:
    return [
        line
        for line in path.read_text().splitlines()
        if line and not line.startswith("#")
    ]


def _write_raw_vcf(path: Path, rows: list[str]) -> None:
    path.write_text(
        "##fileformat=VCFv4.2\n"
        "##contig=<ID=chr2,length=1000000>\n"
        "##contig=<ID=chr1,length=1000000>\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        + "".join(rows)
    )


def _raw_del(chrom: str, pos: int, record_id: str) -> str:
    end = pos + 50
    return (
        f"{chrom}\t{pos}\t{record_id}\tN\t<DEL>\t60\tPASS\t"
        f"SVTYPE=DEL;END={end};SVLEN=-50\n"
    )


def test_correct_output_follows_header_contig_order_then_position(tmp_path):
    raw = tmp_path / "unsorted.vcf"
    out = tmp_path / "sorted.svcf"
    _write_raw_vcf(
        raw,
        [
            _raw_del("chr1", 300, "chr1_300"),
            _raw_del("chr2", 500, "chr2_500"),
            _raw_del("chr2", 100, "chr2_100"),
            _raw_del("chr1", 100, "chr1_100"),
        ],
    )

    result = RUNNER.invoke(
        app,
        ["correct", "-i", str(raw), "-o", str(out)],
    )

    assert result.exit_code == 0, result.output
    ids = [row.split("\t")[2] for row in _data_rows(out)]
    assert ids == ["chr2_100", "chr2_500", "chr1_100", "chr1_300"]


def test_correct_same_position_sort_is_content_deterministic(tmp_path):
    first_raw = tmp_path / "ties_first.vcf"
    second_raw = tmp_path / "ties_second.vcf"
    first_out = tmp_path / "ties_first.svcf"
    second_out = tmp_path / "ties_second.svcf"

    rows = [
        _raw_del("chr2", 100, "zeta"),
        _raw_del("chr2", 100, "alpha"),
    ]
    _write_raw_vcf(first_raw, rows)
    _write_raw_vcf(second_raw, list(reversed(rows)))

    result1 = RUNNER.invoke(
        app,
        ["correct", "-i", str(first_raw), "-o", str(first_out)],
    )
    result2 = RUNNER.invoke(
        app,
        ["correct", "-i", str(second_raw), "-o", str(second_out)],
    )

    assert result1.exit_code == 0, result1.output
    assert result2.exit_code == 0, result2.output
    ids1 = [row.split("\t")[2] for row in _data_rows(first_out)]
    ids2 = [row.split("\t")[2] for row in _data_rows(second_out)]
    assert ids1 == ids2 == ["alpha", "zeta"]


def _sv_event(chrom: str, pos: int, record_id: str) -> SVEvent:
    end = pos + 50
    event = SVEvent(
        chrom,
        pos,
        record_id,
        "N",
        "<DEL>",
        "60",
        "PASS",
        (
            f"SVTYPE=DEL;END={end};SVLEN=50;CHR2={chrom};SUPPORT=5;"
            "SVMETHOD=OctopuSV;RTID=.;AF=.;STRAND=.;RNAMES=."
        ),
        format="GT",
        sample="0/1",
    )
    event.source = "test"
    return event


def test_correct_writer_sort_only_reorders_record_bytes(tmp_path):
    events = [
        _sv_event("chr1", 300, "a"),
        _sv_event("chr2", 500, "b"),
        _sv_event("chr2", 100, "c"),
    ]
    expected_rows = Counter(str(event) for event in events)
    output = tmp_path / "out.svcf"

    write_sv_vcf(
        [
            "##contig=<ID=chr2,length=1000000>",
            "##contig=<ID=chr1,length=1000000>",
        ],
        events,
        output,
    )

    assert Counter(_data_rows(output)) == expected_rows


def _caller_event(chrom: str, pos: int, record_id: str, sources: tuple[str, ...]):
    end = pos + 50
    merged_records = []
    for index, source in enumerate(sources, start=1):
        merged_records.append(
            (
                source,
                "SAMPLE",
                CALLER_FORMAT,
                {
                    "GT": "0/1",
                    "ID": f"{record_id}.src{index}",
                    "SC": f"caller{index}",
                },
            )
        )

    return SimpleNamespace(
        chrom=chrom,
        pos=pos,
        sv_id=record_id,
        ref="N",
        alt="<DEL>",
        quality="60",
        filter="PASS",
        info={
            "SVTYPE": "DEL",
            "END": str(end),
            "SVLEN": "50",
            "CHR2": chrom,
            "SUPPORT": "5",
            "SVMETHOD": "OctopuSV",
            "RTID": ".",
            "AF": ".",
            "STRAND": ".",
            "RNAMES": ".",
        },
        format=CALLER_FORMAT,
        merged_sample_records=merged_records,
        merged_samples=[record[3] for record in merged_records],
        source_file=",".join(sources),
    )


def test_caller_merge_writer_sorts_records_without_reordering_evidence(monkeypatch, tmp_path):
    source_a = str(tmp_path / "A.svcf")
    source_b = str(tmp_path / "B.svcf")
    merger = SVMerger({}, [source_a, source_b])
    monkeypatch.setattr(
        merger,
        "_write_vcf_header",
        lambda handle, contigs, input_files: None,
    )

    two_source = _caller_event("chr2", 100, "chr2_100", (source_b, source_a))
    events = [
        _caller_event("chr1", 300, "chr1_300", (source_a,)),
        _caller_event("chr2", 500, "chr2_500", (source_a,)),
        two_source,
    ]
    output = tmp_path / "caller.svcf"

    merger.write_results(
        output,
        events,
        {"chr2": "1000000", "chr1": "1000000"},
        mode="caller",
        input_files=[source_a, source_b],
    )

    rows = _data_rows(output)
    assert [row.split("\t")[2] for row in rows] == [
        "chr2_100",
        "chr2_500",
        "chr1_300",
    ]

    first = rows[0].split("\t")
    info = dict(
        item.split("=", 1)
        for item in first[7].split(";")
        if "=" in item
    )
    # Evidence ordering remains tied to merge-input order A, B even though the
    # event itself moved to a different outer record position.
    assert info["SOURCES"].split(",") == ["A", "B"]
    assert info["SOURCE_IDS"].split(",") == ["chr2_100.src2", "chr2_100.src1"]
    sample_ids = [sample.split(":")[1] for sample in first[9:]]
    assert sample_ids == ["chr2_100.src2", "chr2_100.src1"]


def _sample_event(chrom: str, pos: int, record_id: str):
    end = pos + 50
    return SimpleNamespace(
        chrom=chrom,
        pos=pos,
        sv_id=record_id,
        ref="N",
        alt="<DEL>",
        quality="60",
        filter="PASS",
        info={
            "SVTYPE": "DEL",
            "END": str(end),
            "SVLEN": "50",
            "CHR2": chrom,
            "SUPPORT": "5",
            "SVMETHOD": "OctopuSV",
            "RTID": ".",
            "AF": ".",
            "STRAND": ".",
            "RNAMES": ".",
        },
        ordered_samples=[
            {
                "GT": "0/1",
                "AD": ".,.",
                "UC": "1",
                "UV": "1",
                "LN": "50",
                "ST": ".",
                "QV": "60",
                "TY": "DEL",
                "ID": f"{record_id}.A",
                "SC": "OctopuSV",
                "REF": "N",
                "ALT": "<DEL>",
                "CO": f"{chrom}_{pos}-{chrom}_{end}",
            },
            {
                "GT": "1/1",
                "AD": ".,.",
                "UC": "1",
                "UV": "1",
                "LN": "50",
                "ST": ".",
                "QV": "60",
                "TY": "DEL",
                "ID": f"{record_id}.B",
                "SC": "OctopuSV",
                "REF": "N",
                "ALT": "<DEL>",
                "CO": f"{chrom}_{pos}-{chrom}_{end}",
            },
        ],
    )


def test_sample_writer_sorts_records_and_keeps_sample_columns_bound(tmp_path):
    mapper = NameMapper(
        [str(tmp_path / "A.svcf"), str(tmp_path / "B.svcf")],
        mode="sample",
        custom_names=["A", "B"],
    )
    writer = MultiSampleWriter(mapper)
    output = tmp_path / "sample.svcf"

    writer.write_results(
        output,
        [
            _sample_event("chr1", 300, "chr1_300"),
            _sample_event("chr2", 500, "chr2_500"),
            _sample_event("chr2", 100, "chr2_100"),
        ],
        {"chr2": "1000000", "chr1": "1000000"},
    )

    rows = _data_rows(output)
    assert [row.split("\t")[2] for row in rows] == [
        "chr2_100",
        "chr2_500",
        "chr1_300",
    ]

    first = rows[0].split("\t")
    format_keys = first[8].split(":")
    id_index = format_keys.index("ID")
    assert first[9].split(":")[id_index] == "chr2_100.A"
    assert first[10].split(":")[id_index] == "chr2_100.B"
