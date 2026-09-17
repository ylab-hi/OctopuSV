from __future__ import annotations

import gzip
import os
from pathlib import Path
from types import SimpleNamespace

import pytest

from octopusv.cli.merge import (
    _preflight_merge_inputs,
    _read_merge_input_list,
)
from octopusv.filtering.svcf_filter import FilterConfig, SVCFFilter
from octopusv.formatter.svcf_to_bed_converter import SVCFtoBEDConverter
from octopusv.formatter.svcf_to_bedpe_converter import SVCFtoBEDPEConverter
from octopusv.formatter.svcf_to_vcf_converter import SVCFtoVCFConverter
from octopusv.merger.name_mapper import NameMapper, default_label_from_path
from octopusv.normalization.contig_normalizer import SVCFContigNormalizer
from octopusv.querying.svcf_query import QueryConfig, SVCFQuery, parse_region
from octopusv.subset.svcf_subset import SVCFSubset, SubsetConfig
from octopusv.utils.atomic_write import atomic_output_path
from octopusv.utils.header_reader import HeaderReader
from octopusv.utils.normal_vcf_parser import parse_vcf
from octopusv.utils.svcf_parser import SVCFEvent, SVCFFileEventCreator
from octopusv.utils.svcf_schema import CALLER_FORMAT, SAMPLE_FORMAT, SVCF_VERSION
from octopusv.utils.svcf_validator import SVCFValidator
from octopusv.utils.text_io import open_text_auto


CALLER_BLOCK = (
    "0/1:5,7:50:+-:60:DEL:caller.1:caller:N:<DEL>:"
    "chr1_100-chr1_150"
)
BASE_INFO = (
    "SVTYPE=DEL;END=150;SVLEN=50;CHR2=chr1;SUPPORT=5;"
    "SVMETHOD=OctopuSV;RTID=.;AF=.;STRAND=+-;RNAMES=."
)


def _caller_svcf_text(
    *,
    info_suffix: str = "",
    record_id: str = "event1",
    contig_length: int = 1000,
    source: str = "caller",
    source_id: str = "caller.1",
) -> str:
    block = CALLER_BLOCK.replace(":caller.1:caller:", f":{source_id}:{source}:")
    info = BASE_INFO + (";" + info_suffix if info_suffix else "")
    return (
        "##fileformat=VCFv4.2\n"
        f"##SVCFVersion={SVCF_VERSION}\n"
        "##OctopuSV_mode=caller\n"
        f"##contig=<ID=chr1,length={contig_length}>\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n"
        f"chr1\t100\t{record_id}\tN\t<DEL>\t60\tPASS\t{info}\t"
        f"{CALLER_FORMAT}\t{block}\n"
    )


def _merged_caller_svcf_text(*, second_record: str | None = None) -> str:
    text = _caller_svcf_text(
        info_suffix="SOURCES=caller;SOURCE_IDS=caller.1",
        source="caller",
        source_id="caller.1",
    )
    if second_record is not None:
        text += second_record
    return text


def _write_gzip(path: Path, text: str) -> None:
    with gzip.open(path, "wt", encoding="utf-8") as handle:
        handle.write(text)


def _write_multimember_gzip(path: Path, parts: list[str]) -> None:
    payload = b""
    for part in parts:
        payload += gzip.compress(part.encode("utf-8"))
    path.write_bytes(payload)


@pytest.mark.parametrize(
    "kind",
    ["plain", "gzip_suffix", "gzip_no_suffix", "plain_named_gz", "bom", "multi_member", "empty"],
)
def test_open_text_auto_content_detection(kind, tmp_path):
    path = tmp_path / ("input.gz" if kind in {"gzip_suffix", "plain_named_gz"} else "input")
    expected = "hello\n"

    if kind == "plain":
        path.write_text(expected)
    elif kind == "gzip_suffix":
        _write_gzip(path, expected)
    elif kind == "gzip_no_suffix":
        _write_gzip(path, expected)
    elif kind == "plain_named_gz":
        path.write_text(expected)
    elif kind == "bom":
        path.write_bytes(b"\xef\xbb\xbfhello\n")
    elif kind == "multi_member":
        _write_multimember_gzip(path, ["hel", "lo\n"])
    else:
        path.write_bytes(b"")
        expected = ""

    with open_text_auto(path) as handle:
        assert handle.read() == expected


def test_correct_parser_accepts_gzip_without_gzip_suffix(tmp_path):
    raw = (
        "##fileformat=VCFv4.2\n"
        "##source=TestCaller\n"
        "##contig=<ID=chr1,length=1000>\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n"
        "chr1\t100\tr1\tN\t<DEL>\t60\tPASS\tSVTYPE=DEL;END=150;SVLEN=-50\tGT:AD\t0/1:5,7\n"
    )
    path = tmp_path / "input.vcf"
    _write_gzip(path, raw)

    contigs, same_chr, diff_chr, non_bnd = parse_vcf(path)
    assert contigs == ["##contig=<ID=chr1,length=1000>"]
    assert same_chr == []
    assert diff_chr == []
    assert len(non_bnd) == 1
    assert non_bnd[0].id == "r1"


@pytest.mark.parametrize(
    ("path", "expected"),
    [
        ("/x/x.svcf", "x"),
        ("/x/x.svcf.gz", "x"),
        ("/x/x.vcf.bgz", "x"),
        ("/x/a.v2.svcf", "a.v2"),
        ("/x/X.SVCF.GZ", "X"),
    ],
)
def test_default_label_from_path_table(path, expected):
    assert default_label_from_path(path) == expected


def test_name_mapper_custom_names_never_fall_back_to_basename():
    mapper = NameMapper(["/data/a.svcf"], custom_names=["A"])
    with pytest.raises(ValueError, match="refusing to fall back"):
        mapper.get_display_name("/different/a.svcf")


def test_input_list_relative_paths_preserve_symlink_identity(tmp_path):
    store = tmp_path / "store"
    store.mkdir()
    target = store / "abc123.svcf.gz"
    _write_gzip(target, _caller_svcf_text())

    lists = tmp_path / "lists"
    lists.mkdir()
    visible = lists / "patient007.svcf.gz"
    visible.symlink_to(target)
    input_list = lists / "cohort.tsv"
    input_list.write_text("patient007.svcf.gz\n")

    paths, labels = _read_merge_input_list(input_list)
    assert labels is None
    assert paths == [visible.absolute()]
    assert paths[0].name == "patient007.svcf.gz"
    assert paths[0].resolve().name == "abc123.svcf.gz"
    assert default_label_from_path(paths[0]) == "patient007"


@pytest.mark.parametrize(
    "content,match",
    [
        ("a.svcf\tA\textra\n", "PATH or PATH<TAB>LABEL"),
        ("a.svcf\tA\nb.svcf\n", "must provide a TAB-separated label"),
        ("a.svcf\t   \n", "label is empty"),
    ],
)
def test_input_list_rejects_ambiguous_rows(tmp_path, content, match):
    (tmp_path / "a.svcf").write_text(_caller_svcf_text())
    (tmp_path / "b.svcf").write_text(_caller_svcf_text(record_id="b"))
    list_path = tmp_path / "bad.tsv"
    list_path.write_text(content)

    with pytest.raises(ValueError, match=match):
        _read_merge_input_list(list_path)


def test_input_list_supports_spaces_comments_crlf_and_gzip(tmp_path):
    data = tmp_path / "data dir"
    data.mkdir()
    a = data / "sample one.svcf"
    a.write_text(_caller_svcf_text())
    list_path = tmp_path / "inputs.tsv"
    _write_gzip(list_path, "  # comment\r\n\r\ndata dir/sample one.svcf\tSampleOne\r\n")

    paths, labels = _read_merge_input_list(list_path)
    assert paths == [a.absolute()]
    assert labels == ["SampleOne"]


def test_invalid_label_fails_before_shape_scan(tmp_path, monkeypatch):
    path = tmp_path / "patient 1.svcf"
    path.write_text("this is intentionally not valid SVCF\n")

    import octopusv.cli.merge as merge_module

    def should_not_run(*args, **kwargs):
        raise AssertionError("shape scan must not run before label validation")

    monkeypatch.setattr(merge_module, "_preflight_svcf_shape", should_not_run)

    with pytest.raises(ValueError, match="Invalid caller label"):
        _preflight_merge_inputs(
            input_files=[path],
            labels=["patient 1"],
            mode="caller",
        )


@pytest.mark.parametrize("label", ["bad name", "bad,name", "bad;name", "bad=name", ".", ""])
def test_explicit_label_reserved_character_matrix(tmp_path, monkeypatch, label):
    path = tmp_path / "x.svcf"
    path.write_text("invalid on purpose\n")

    import octopusv.cli.merge as merge_module

    monkeypatch.setattr(
        merge_module,
        "_preflight_svcf_shape",
        lambda *args, **kwargs: (_ for _ in ()).throw(
            AssertionError("shape scan must not run before label validation")
        ),
    )

    with pytest.raises(ValueError, match="Invalid caller label"):
        _preflight_merge_inputs(input_files=[path], labels=[label], mode="caller")


def test_preflight_rejects_conflicting_contig_lengths(tmp_path):
    a = tmp_path / "a.svcf"
    b = tmp_path / "b.svcf"
    a.write_text(_caller_svcf_text(contig_length=1000))
    b.write_text(_caller_svcf_text(contig_length=2000, record_id="b"))

    with pytest.raises(ValueError, match="Conflicting ##contig lengths"):
        _preflight_merge_inputs(
            input_files=[a, b],
            labels=["a", "b"],
            mode="sample",
        )


def test_parser_rejects_malformed_record_with_line_number(tmp_path):
    path = tmp_path / "bad.svcf"
    path.write_text(
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n"
        "chr1\t100\ttoo_short\n"
    )

    creator = SVCFFileEventCreator([str(path)])
    with pytest.raises(ValueError, match=r"line 3"):
        creator.parse()


def test_co_fallback_uses_shared_parser_for_hyphenated_contigs():
    sample = (
        "0/1:5,7:.:.:60:DEL:e1:caller:N:<DEL>:"
        "chr1-alt_100-chr2-alt_200"
    )
    event = SVCFEvent(
        "chr1-alt",
        "100",
        "e1",
        "N",
        "<DEL>",
        "60",
        "PASS",
        "SVTYPE=DEL;END=.;SVLEN=.;CHR2=chr2-alt;SUPPORT=5",
        CALLER_FORMAT,
        sample,
        source_file="x.svcf",
        sample_name="SAMPLE",
    )
    assert (event.start_chrom, event.start_pos, event.end_chrom, event.end_pos) == (
        "chr1-alt",
        100,
        "chr2-alt",
        200,
    )


def test_bed_bnd_is_one_bp_and_keeps_bnd_type_name():
    event = SimpleNamespace(
        chrom="chr1",
        pos=100,
        sv_id="bnd1",
        sv_type="BND",
        info={"CHR2": "chr2", "END": "5000000", "SVLEN": ".", "STRAND": "+-"},
        quality="60",
    )
    fields = SVCFtoBEDConverter([event]).convert().splitlines()[-1].split("\t")
    assert fields[:3] == ["chr1", "99", "100"]
    assert "_BND_" in fields[3]
    assert "_TRA_" not in fields[3]


def test_bedpe_paired_bnd_uses_explicit_mate():
    event = SimpleNamespace(
        chrom="chr1",
        pos=100,
        sv_id="bnd1",
        sv_type="BND",
        info={"CHR2": "chr2", "END": "250", "SVLEN": ".", "STRAND": "+-"},
        quality="60",
    )
    fields = SVCFtoBEDPEConverter([event]).convert().splitlines()[-1].split("\t")
    assert fields[:6] == ["chr1", "99", "100", "chr2", "249", "250"]


def test_bedpe_single_breakend_fails_instead_of_inventing_mate():
    event = SimpleNamespace(
        chrom="chr1",
        pos=100,
        sv_id="single",
        sv_type="BND",
        info={"CHR2": ".", "END": ".", "SVLEN": ".", "STRAND": "."},
        quality="60",
    )
    with pytest.raises(ValueError, match="requires explicit CHR2"):
        SVCFtoBEDPEConverter([event]).convert()


def test_filter_source_prefers_sources_over_misleading_id(tmp_path):
    input_path = tmp_path / "merged.svcf"
    output_path = tmp_path / "out.svcf"
    input_path.write_text(
        _caller_svcf_text(
            record_id="pbsv.DEL.1",
            info_suffix="SOURCES=sniffles;SOURCE_IDS=caller.1",
            source="cuteSV",
        )
    )

    kept = SVCFFilter(
        FilterConfig(
            input_file=str(input_path),
            output_file=str(output_path),
            sources={"sniffles"},
        )
    ).run()
    assert kept["output_records"] == 1
    assert kept["source_detection_used"] == {"INFO/SOURCES": 1}

    dropped = SVCFFilter(
        FilterConfig(
            input_file=str(input_path),
            output_file=str(tmp_path / "drop.svcf"),
            sources={"pbsv"},
        )
    ).run()
    assert dropped["output_records"] == 0


def test_filter_single_evidence_uses_explicit_sc_not_id_prefix(tmp_path):
    input_path = tmp_path / "single.svcf"
    input_path.write_text(
        _caller_svcf_text(record_id="pbsv.DEL.1", source="cuteSV")
    )
    summary = SVCFFilter(
        FilterConfig(
            input_file=str(input_path),
            output_file=str(tmp_path / "out.svcf"),
            sources={"cuteSV"},
        )
    ).run()
    assert summary["output_records"] == 1
    assert summary["source_detection_used"] == {"FORMAT/SC": 1}


def test_filter_source_fails_when_no_explicit_source_exists(tmp_path):
    path = tmp_path / "ambiguous.svcf"
    path.write_text(
        _caller_svcf_text(record_id="pbsv.DEL.1", source=".")
    )
    with pytest.raises(ValueError, match="does not infer source identity"):
        SVCFFilter(
            FilterConfig(
                input_file=str(path),
                output_file=str(tmp_path / "out.svcf"),
                sources={"pbsv"},
            )
        ).run()


def test_atomic_output_path_preserves_old_file_and_permissions(tmp_path):
    output = tmp_path / "out.txt"
    output.write_text("old\n")
    output.chmod(0o640)

    with pytest.raises(RuntimeError, match="boom"):
        with atomic_output_path(output) as temp:
            temp.write_text("partial\n")
            raise RuntimeError("boom")

    assert output.read_text() == "old\n"
    assert output.stat().st_mode & 0o777 == 0o640
    assert not list(tmp_path.glob(".out.txt.*.tmp"))

    with atomic_output_path(output) as temp:
        temp.write_text("new\n")

    assert output.read_text() == "new\n"
    assert output.stat().st_mode & 0o777 == 0o640



def test_atomic_output_path_can_replace_read_only_target_and_preserve_mode(tmp_path):
    output = tmp_path / "readonly.txt"
    output.write_text("old\n")
    output.chmod(0o444)

    with atomic_output_path(output) as temp:
        Path(temp).write_text("new\n")

    assert output.read_text() == "new\n"
    assert output.stat().st_mode & 0o777 == 0o444

def _truncated_svcf_text() -> str:
    return _merged_caller_svcf_text(second_record="chr1\t300\ttruncated\n")


def _run_filter(input_path: Path, output_path: Path):
    return SVCFFilter(FilterConfig(input_file=str(input_path), output_file=str(output_path))).run()


def _run_subset(input_path: Path, output_path: Path):
    return SVCFSubset(
        SubsetConfig(
            input_file=str(input_path),
            output_file=str(output_path),
            mode="caller",
            selected_callers=["caller"],
        )
    ).run()


def _run_query(input_path: Path, output_path: Path):
    return SVCFQuery(
        QueryConfig(
            input_file=str(input_path),
            output_file=str(output_path),
            targets=[parse_region("chr1:1-1000")],
        )
    ).run()


def _run_normalize(input_path: Path, output_path: Path):
    return SVCFContigNormalizer(input_path, output_path, "chr").run()


def _run_svcf2vcf(input_path: Path, output_path: Path):
    return SVCFtoVCFConverter(events=None, input_svcf_file=str(input_path)).convert_to_file(output_path)


@pytest.mark.parametrize(
    "runner",
    [_run_filter, _run_subset, _run_query, _run_normalize, _run_svcf2vcf],
    ids=["filter", "subset", "query", "normalize", "svcf2vcf"],
)
def test_streaming_writers_fail_loud_and_leave_existing_output_unchanged(tmp_path, runner):
    input_path = tmp_path / "truncated.svcf"
    input_path.write_text(_truncated_svcf_text())
    output = tmp_path / "out.txt"
    output.write_text("ORIGINAL\n")

    with pytest.raises(ValueError, match=r"line 7"):
        runner(input_path, output)

    assert output.read_text() == "ORIGINAL\n"
    assert not list(tmp_path.glob(f".{output.name}.*.tmp"))


def test_header_reader_keeps_legacy_marker_key_semantics(tmp_path):
    caller = tmp_path / "caller.svcf"
    caller.write_text(_caller_svcf_text())
    reader = HeaderReader(caller)
    reader.read()
    contract = reader.to_contract()
    assert contract["declared_mode"] == "caller"
    assert contract["has_octopusv_mode_marker"] is False

    multi = tmp_path / "legacy_multi.svcf"
    multi.write_text(
        "##fileformat=VCFv4.2\n"
        "##OctopuSV_mode=multi\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\n"
    )
    reader = HeaderReader(multi)
    reader.read()
    assert reader.to_contract()["has_octopusv_mode_marker"] is True


def test_validator_reads_gzip_without_gzip_suffix(tmp_path):
    path = tmp_path / "caller.data"
    _write_gzip(path, _caller_svcf_text())
    validator = SVCFValidator(str(path))
    validator.validate()
    assert validator.errors == []


def test_spec_schema_constants_do_not_drift():
    spec = Path(__file__).parents[1] / "docs" / "SVCF_specifications.md"
    text = spec.read_text(encoding="utf-8")
    assert f"SVCFVersion={SVCF_VERSION}" in text
    assert CALLER_FORMAT in text
    assert SAMPLE_FORMAT in text
