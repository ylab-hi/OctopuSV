from io import StringIO
from pathlib import Path
from types import SimpleNamespace

import pytest

from octopusv.cli.merge import _preflight_merge_inputs
from octopusv.merger.multi_sample_writer import MultiSampleWriter
from octopusv.merger.name_mapper import NameMapper
from octopusv.merger.sv_merge_writer import MergeWriterMixin
from octopusv.utils.svcf_schema import (
    CALLER_FORMAT,
    SAMPLE_FORMAT,
    SVCF_VERSION,
    validate_source_id,
    validate_source_label,
)
from octopusv.utils.svcf_validator import SVCFValidator


CALLER_BLOCK = (
    "0/1:5,7:10:.:60:INS:caller.1:caller:N:<INS>:"
    "chr1_100-chr1_110"
)
SAMPLE_BLOCK = (
    "0/1:.,.:1:1:10:.:60:INS:sample.1:OctopuSV:N:<INS>:"
    "chr1_100-chr1_110"
)
BASE_INFO = (
    "SVTYPE=INS;END=110;SVLEN=10;CHR2=chr1;SUPPORT=5;"
    "SVMETHOD=OctopuSV;RTID=.;AF=.;STRAND=.;RNAMES=."
)


def _write_record_file(
    path: Path,
    *,
    meta: list[str],
    fmt: str,
    blocks: list[str],
    samples: list[str] | None = None,
    info_suffix: str = "",
) -> Path:
    if samples is None:
        samples = ["SAMPLE"]

    header = (
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"
        + "\t".join(samples)
    )
    info = BASE_INFO + (";" + info_suffix if info_suffix else "")
    record = (
        f"chr1\t100\te1\tN\t<INS>\t60\tPASS\t{info}\t{fmt}\t"
        + "\t".join(blocks)
    )
    path.write_text("\n".join(meta + [header, record]) + "\n")
    return path


class _Writer(MergeWriterMixin):
    pass


def test_caller_merge_writer_declares_svcf_11_caller_mode(tmp_path):
    input_file = tmp_path / "input.svcf"
    input_file.write_text(
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n"
    )

    out = StringIO()
    _Writer()._write_vcf_header(out, {"chr1": 1000}, [str(input_file)])
    text = out.getvalue()

    assert f"##SVCFVersion={SVCF_VERSION}\n" in text
    assert "##OctopuSV_mode=caller\n" in text


def test_sample_writer_declares_svcf_11_multi_mode(tmp_path):
    mapper = NameMapper([str(tmp_path / "s1.svcf")], mode="sample")
    writer = MultiSampleWriter(mapper)
    out = StringIO()

    writer._write_header(out, {"chr1": 1000})
    text = out.getvalue()

    assert f"##SVCFVersion={SVCF_VERSION}\n" in text
    assert "##OctopuSV_mode=multi\n" in text


def test_validator_accepts_versioned_caller_schema(tmp_path):
    path = _write_record_file(
        tmp_path / "caller.svcf",
        meta=[
            "##fileformat=VCFv4.2",
            f"##SVCFVersion={SVCF_VERSION}",
            "##OctopuSV_mode=caller",
        ],
        fmt=CALLER_FORMAT,
        blocks=[CALLER_BLOCK],
    )

    validator = SVCFValidator(str(path))
    validator.validate()

    assert validator.svcf_version == SVCF_VERSION
    assert validator.declared_mode == "caller"
    assert validator.errors == []


def test_validator_accepts_versioned_multi_schema(tmp_path):
    path = _write_record_file(
        tmp_path / "multi.svcf",
        meta=[
            "##fileformat=VCFv4.2",
            f"##SVCFVersion={SVCF_VERSION}",
            "##OctopuSV_mode=multi",
        ],
        fmt=SAMPLE_FORMAT,
        blocks=[SAMPLE_BLOCK],
        samples=["sample1"],
        info_suffix="SOURCES=sample1;SOURCE_IDS=sample.1",
    )

    validator = SVCFValidator(str(path))
    validator.validate()

    assert validator.svcf_version == SVCF_VERSION
    assert validator.declared_mode == "multi"
    assert validator.errors == []


def test_versioned_file_requires_explicit_mode(tmp_path):
    path = _write_record_file(
        tmp_path / "missing_mode.svcf",
        meta=["##fileformat=VCFv4.2", f"##SVCFVersion={SVCF_VERSION}"],
        fmt=CALLER_FORMAT,
        blocks=[CALLER_BLOCK],
    )

    validator = SVCFValidator(str(path))
    validator.validate()

    assert "E_VER_001" in {issue.code for issue in validator.errors}


def test_versioned_caller_marker_rejects_sample_schema(tmp_path):
    path = _write_record_file(
        tmp_path / "bad_schema.svcf",
        meta=[
            "##fileformat=VCFv4.2",
            f"##SVCFVersion={SVCF_VERSION}",
            "##OctopuSV_mode=caller",
        ],
        fmt=SAMPLE_FORMAT,
        blocks=[SAMPLE_BLOCK],
    )

    validator = SVCFValidator(str(path))
    validator.validate()

    assert "E_FMT_001" in {issue.code for issue in validator.errors}


def test_merge_preflight_rejects_ordinary_vcf_before_parsing(tmp_path):
    path = _write_record_file(
        tmp_path / "ordinary.vcf",
        meta=["##fileformat=VCFv4.2"],
        fmt="GT:AD:DP",
        blocks=["0/1:5,7:12"],
    )

    with pytest.raises(ValueError, match="not recognized as SVCF"):
        _preflight_merge_inputs(
            input_files=[path],
            labels=["ordinary"],
            mode="caller",
        )


def test_merge_preflight_rejects_sample_schema_even_if_marker_was_stripped(tmp_path):
    path = _write_record_file(
        tmp_path / "stripped_multi.svcf",
        meta=["##fileformat=VCFv4.2"],
        fmt=SAMPLE_FORMAT,
        blocks=[SAMPLE_BLOCK],
        samples=["sample1"],
    )

    with pytest.raises(ValueError, match="sample-mode SVCF schema"):
        _preflight_merge_inputs(
            input_files=[path],
            labels=["sample1"],
            mode="sample",
        )


def test_sample_merge_accepts_versioned_caller_merged_input(tmp_path):
    path = _write_record_file(
        tmp_path / "caller_merged.svcf",
        meta=[
            "##fileformat=VCFv4.2",
            f"##SVCFVersion={SVCF_VERSION}",
            "##OctopuSV_mode=caller",
        ],
        fmt=CALLER_FORMAT,
        blocks=[CALLER_BLOCK, CALLER_BLOCK],
        info_suffix="SOURCES=a,b;SOURCE_IDS=caller.1,caller.1",
    )

    _preflight_merge_inputs(
        input_files=[path],
        labels=["sample1"],
        mode="sample",
    )


def _write_two_sample_v11(path: Path) -> Path:
    return _write_record_file(
        path,
        meta=[
            "##fileformat=VCFv4.2",
            f"##SVCFVersion={SVCF_VERSION}",
            "##OctopuSV_mode=multi",
        ],
        fmt=SAMPLE_FORMAT,
        blocks=[SAMPLE_BLOCK, SAMPLE_BLOCK],
        samples=["sample1", "sample2"],
        info_suffix="SOURCES=sample1,sample2;SOURCE_IDS=sample.1,sample.1",
    )


def _assert_v11_identity_preserved(path: Path):
    text = path.read_text()
    assert f"##SVCFVersion={SVCF_VERSION}\n" in text
    assert "##OctopuSV_mode=multi\n" in text

    validator = SVCFValidator(str(path))
    validator.validate()
    assert validator.errors == []


def test_filter_preserves_svcf_11_identity_headers(tmp_path):
    from octopusv.filtering.svcf_filter import FilterConfig, SVCFFilter

    input_path = _write_two_sample_v11(tmp_path / "input.svcf")
    output_path = tmp_path / "filtered.svcf"

    SVCFFilter(
        FilterConfig(
            input_file=str(input_path),
            output_file=str(output_path),
        )
    ).run()

    _assert_v11_identity_preserved(output_path)


def test_normalize_contigs_preserves_svcf_11_identity_headers(tmp_path):
    from octopusv.normalization.contig_normalizer import SVCFContigNormalizer

    input_path = _write_two_sample_v11(tmp_path / "input.svcf")
    output_path = tmp_path / "normalized.svcf"

    SVCFContigNormalizer(
        input_file=input_path,
        output_file=output_path,
        style="chr",
    ).run()

    _assert_v11_identity_preserved(output_path)


def test_subset_preserves_svcf_11_identity_headers(tmp_path):
    from octopusv.subset.svcf_subset import SubsetConfig, SVCFSubset

    input_path = _write_two_sample_v11(tmp_path / "input.svcf")
    output_path = tmp_path / "subset.svcf"

    SVCFSubset(
        SubsetConfig(
            input_file=str(input_path),
            output_file=str(output_path),
            selected_samples=["sample1"],
        )
    ).run()

    _assert_v11_identity_preserved(output_path)


def test_svcf2vcf_rejects_unsupported_version(tmp_path):
    from octopusv.formatter.svcf_to_vcf_converter import SVCFtoVCFConverter

    path = _write_record_file(
        tmp_path / "future.svcf",
        meta=[
            "##fileformat=VCFv4.2",
            "##SVCFVersion=1.2",
            "##OctopuSV_mode=caller",
        ],
        fmt=CALLER_FORMAT,
        blocks=[CALLER_BLOCK],
    )

    with pytest.raises(ValueError, match="Unsupported SVCFVersion"):
        SVCFtoVCFConverter(input_svcf_file=path)


def test_svcf2vcf_rejects_versioned_mode_schema_conflict(tmp_path):
    from octopusv.formatter.svcf_to_vcf_converter import SVCFtoVCFConverter

    path = _write_record_file(
        tmp_path / "bad_consumer_schema.svcf",
        meta=[
            "##fileformat=VCFv4.2",
            f"##SVCFVersion={SVCF_VERSION}",
            "##OctopuSV_mode=caller",
        ],
        fmt=SAMPLE_FORMAT,
        blocks=[SAMPLE_BLOCK],
    )

    with pytest.raises(ValueError, match="requires FORMAT"):
        SVCFtoVCFConverter(input_svcf_file=path)


def test_stat_reader_rejects_unsupported_version(tmp_path):
    from octopusv.stater.stat_reader import read_records

    path = _write_record_file(
        tmp_path / "future_stat.svcf",
        meta=[
            "##fileformat=VCFv4.2",
            "##SVCFVersion=1.2",
            "##OctopuSV_mode=caller",
        ],
        fmt=CALLER_FORMAT,
        blocks=[CALLER_BLOCK],
    )

    with pytest.raises(ValueError, match="Unsupported SVCFVersion"):
        read_records(path)


def test_stat_reader_rejects_versioned_mode_schema_conflict(tmp_path):
    from octopusv.stater.stat_reader import read_records

    path = _write_record_file(
        tmp_path / "bad_stat_schema.svcf",
        meta=[
            "##fileformat=VCFv4.2",
            f"##SVCFVersion={SVCF_VERSION}",
            "##OctopuSV_mode=multi",
        ],
        fmt=CALLER_FORMAT,
        blocks=[CALLER_BLOCK],
        samples=["sample1"],
    )

    with pytest.raises(ValueError, match="requires FORMAT"):
        read_records(path)


def test_sample_merge_preflight_checks_every_record_format(tmp_path):
    path = tmp_path / "mixed_records.svcf"
    header = (
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE"
    )
    caller_record = (
        f"chr1\t100\te1\tN\t<INS>\t60\tPASS\t{BASE_INFO}\t"
        f"{CALLER_FORMAT}\t{CALLER_BLOCK}"
    )
    sample_record = (
        f"chr1\t200\te2\tN\t<INS>\t60\tPASS\t{BASE_INFO}\t"
        f"{SAMPLE_FORMAT}\t{SAMPLE_BLOCK}"
    )
    path.write_text(
        "##fileformat=VCFv4.2\n" + header + "\n" + caller_record + "\n" + sample_record + "\n"
    )

    with pytest.raises(ValueError, match="sample-mode SVCF schema"):
        _preflight_merge_inputs(
            input_files=[path],
            labels=["sample1"],
            mode="sample",
        )


def test_caller_header_generation_fails_loudly_without_fallback(monkeypatch):
    writer = _Writer()
    out = StringIO()

    def boom(_input_files):
        raise RuntimeError("header extraction failed")

    monkeypatch.setattr(writer, "_extract_and_merge_all_headers", boom)

    with pytest.raises(RuntimeError, match="header extraction failed"):
        writer._write_vcf_header(out, {"chr1": 1000}, ["input.svcf"])

    # Extraction happens before any identity/header bytes are written, so a
    # failed writer cannot leave a plausible-but-corrupt partial SVCF header.
    assert out.getvalue() == ""


def test_versioned_sample_writer_preserves_info_flag_syntax(tmp_path):
    mapper = NameMapper([str(tmp_path / "s1.svcf")], mode="sample")
    writer = MultiSampleWriter(mapper)
    event = SimpleNamespace(
        chrom="chr1",
        pos=100,
        sv_id="e1",
        ref="N",
        alt="<INS>",
        quality=60,
        filter="PASS",
        info={
            "SVTYPE": "INS",
            "END": 110,
            "SVLEN": 10,
            "CHR2": "chr1",
            "SUPPORT": 5,
            "SVMETHOD": "OctopuSV",
            "RTID": ".",
            "AF": ".",
            "STRAND": ".",
            "RNAMES": ".",
            "PRECISE": True,
        },
        ordered_samples=[{
            "GT": "0/1",
            "AD": ".,.",
            "UC": "1",
            "UV": "1",
            "LN": "10",
            "ST": ".",
            "QV": "60",
            "TY": "INS",
            "ID": "sample.1",
            "SC": "OctopuSV",
            "REF": "N",
            "ALT": "<INS>",
            "CO": "chr1_100-chr1_110",
        }],
    )

    out = StringIO()
    writer._write_event(out, event)
    info = out.getvalue().rstrip("\n").split("\t")[7]

    assert "PRECISE" in info.split(";")
    assert "PRECISE=True" not in info


def test_caller_header_extraction_does_not_silently_skip_failed_input(monkeypatch):
    def boom(_path):
        raise RuntimeError("cannot read input header")

    monkeypatch.setattr(
        "octopusv.utils.svcf_utils.extract_original_header_definitions",
        boom,
    )

    with pytest.raises(RuntimeError, match="cannot read input header"):
        _Writer()._extract_and_merge_all_headers(["broken.svcf"])


def test_validator_rejects_versioned_caller_with_multiple_header_columns(tmp_path):
    path = _write_record_file(
        tmp_path / "caller_two_header_cols.svcf",
        meta=[
            "##fileformat=VCFv4.2",
            f"##SVCFVersion={SVCF_VERSION}",
            "##OctopuSV_mode=caller",
        ],
        fmt=CALLER_FORMAT,
        blocks=[CALLER_BLOCK, CALLER_BLOCK],
        samples=["sample1", "sample2"],
        info_suffix="SOURCES=a,b;SOURCE_IDS=caller.1,caller.1",
    )

    validator = SVCFValidator(str(path))
    validator.validate()

    assert "E_MODE_001" in {issue.code for issue in validator.errors}


def test_merge_preflight_rejects_versioned_caller_with_multiple_header_columns(tmp_path):
    path = _write_record_file(
        tmp_path / "caller_two_header_cols.svcf",
        meta=[
            "##fileformat=VCFv4.2",
            f"##SVCFVersion={SVCF_VERSION}",
            "##OctopuSV_mode=caller",
        ],
        fmt=CALLER_FORMAT,
        blocks=[CALLER_BLOCK, CALLER_BLOCK],
        samples=["sample1", "sample2"],
        info_suffix="SOURCES=a,b;SOURCE_IDS=caller.1,caller.1",
    )

    with pytest.raises(ValueError, match="requires exactly one #CHROM trailing column"):
        _preflight_merge_inputs(
            input_files=[path],
            labels=["sample1"],
            mode="sample",
        )


def test_svcf2vcf_rejects_versioned_caller_with_multiple_header_columns(tmp_path):
    from octopusv.formatter.svcf_to_vcf_converter import SVCFtoVCFConverter

    path = _write_record_file(
        tmp_path / "caller_two_header_cols.svcf",
        meta=[
            "##fileformat=VCFv4.2",
            f"##SVCFVersion={SVCF_VERSION}",
            "##OctopuSV_mode=caller",
        ],
        fmt=CALLER_FORMAT,
        blocks=[CALLER_BLOCK, CALLER_BLOCK],
        samples=["sample1", "sample2"],
        info_suffix="SOURCES=a,b;SOURCE_IDS=caller.1,caller.1",
    )

    with pytest.raises(ValueError, match="requires exactly one #CHROM trailing column"):
        SVCFtoVCFConverter(input_svcf_file=path)


def test_stat_reader_rejects_versioned_caller_with_multiple_header_columns(tmp_path):
    from octopusv.stater.stat_reader import read_records

    path = _write_record_file(
        tmp_path / "caller_two_header_cols.svcf",
        meta=[
            "##fileformat=VCFv4.2",
            f"##SVCFVersion={SVCF_VERSION}",
            "##OctopuSV_mode=caller",
        ],
        fmt=CALLER_FORMAT,
        blocks=[CALLER_BLOCK, CALLER_BLOCK],
        samples=["sample1", "sample2"],
        info_suffix="SOURCES=a,b;SOURCE_IDS=caller.1,caller.1",
    )

    with pytest.raises(ValueError, match="requires exactly one #CHROM trailing column"):
        read_records(path)


def _write_duplicate_source_caller_v11(path: Path) -> Path:
    block_1 = CALLER_BLOCK
    block_2 = (
        "0/1:6,8:10:.:55:INS:caller.2:caller:N:<INS>:"
        "chr1_100-chr1_110"
    )
    return _write_record_file(
        path,
        meta=[
            "##fileformat=VCFv4.2",
            f"##SVCFVersion={SVCF_VERSION}",
            "##OctopuSV_mode=caller",
        ],
        fmt=CALLER_FORMAT,
        blocks=[block_1, block_2],
        samples=["SAMPLE"],
        info_suffix="SOURCES=caller,caller;SOURCE_IDS=caller.1,caller.2",
    )


def _assert_v11_caller_contract(path: Path):
    text = path.read_text()
    assert f"##SVCFVersion={SVCF_VERSION}\n" in text
    assert "##OctopuSV_mode=caller\n" in text

    validator = SVCFValidator(str(path))
    validator.validate()
    assert validator.errors == []


def test_filter_preserves_versioned_caller_evidence_contract(tmp_path):
    from octopusv.filtering.svcf_filter import FilterConfig, SVCFFilter

    input_path = _write_duplicate_source_caller_v11(tmp_path / "caller.svcf")
    output_path = tmp_path / "caller.filtered.svcf"

    SVCFFilter(
        FilterConfig(
            input_file=str(input_path),
            output_file=str(output_path),
        )
    ).run()

    _assert_v11_caller_contract(output_path)


def test_normalize_contigs_preserves_versioned_caller_evidence_contract(tmp_path):
    from octopusv.normalization.contig_normalizer import SVCFContigNormalizer

    input_path = _write_duplicate_source_caller_v11(tmp_path / "caller.svcf")
    output_path = tmp_path / "caller.normalized.svcf"

    SVCFContigNormalizer(
        input_file=input_path,
        output_file=output_path,
        style="chr",
    ).run()

    _assert_v11_caller_contract(output_path)


def test_subset_preserves_versioned_caller_evidence_contract(tmp_path):
    from octopusv.subset.svcf_subset import SubsetConfig, SVCFSubset

    input_path = _write_duplicate_source_caller_v11(tmp_path / "caller.svcf")
    output_path = tmp_path / "caller.subset.svcf"

    SVCFSubset(
        SubsetConfig(
            input_file=str(input_path),
            output_file=str(output_path),
            selected_callers=["caller"],
        )
    ).run()

    _assert_v11_caller_contract(output_path)


@pytest.mark.parametrize("unsafe", ["bad,name", "bad;name", "bad=name", "bad name"])
def test_svcf11_source_atom_contract_rejects_reserved_characters(unsafe):
    with pytest.raises(ValueError, match="cannot be represented safely"):
        validate_source_label(unsafe)

    with pytest.raises(ValueError, match="cannot be represented safely"):
        validate_source_id(unsafe)


def test_svcf11_source_atom_contract_reserves_dot_for_missing_ids_only():
    assert validate_source_id(".") == "."
    with pytest.raises(ValueError, match="reserved for missing values"):
        validate_source_label(".")


def _minimal_caller_event():
    return SimpleNamespace(
        info={},
        format=CALLER_FORMAT,
    )


@pytest.mark.parametrize("unsafe", ["bad,name", "bad;name", "bad=name", "bad name"])
def test_caller_writer_rejects_unsafe_source_labels(monkeypatch, tmp_path, unsafe):
    writer = _Writer()
    monkeypatch.setattr(writer, "_write_vcf_header", lambda *args, **kwargs: None)
    monkeypatch.setattr(
        writer,
        "_prepare_caller_records",
        lambda *args, **kwargs: [
            {
                "source_name": unsafe,
                "source_id": "caller.1",
                "sample_data": {},
            }
        ],
    )

    with pytest.raises(ValueError, match="SOURCES item"):
        writer.write_results(
            tmp_path / "caller_bad_source.svcf",
            [_minimal_caller_event()],
            {},
            mode="caller",
            input_files=[],
        )


@pytest.mark.parametrize("unsafe", ["bad,id", "bad;id", "bad=id", "bad id"])
def test_caller_writer_rejects_unsafe_source_ids(monkeypatch, tmp_path, unsafe):
    writer = _Writer()
    monkeypatch.setattr(writer, "_write_vcf_header", lambda *args, **kwargs: None)
    monkeypatch.setattr(
        writer,
        "_prepare_caller_records",
        lambda *args, **kwargs: [
            {
                "source_name": "caller",
                "source_id": unsafe,
                "sample_data": {},
            }
        ],
    )

    with pytest.raises(ValueError, match="SOURCE_IDS item"):
        writer.write_results(
            tmp_path / "caller_bad_id.svcf",
            [_minimal_caller_event()],
            {},
            mode="caller",
            input_files=[],
        )


def _sample_mode_event(source_id="sample.1"):
    return SimpleNamespace(
        chrom="chr1",
        pos=100,
        sv_id="event.1",
        ref="N",
        alt="<INS>",
        quality=60,
        filter="PASS",
        info={},
        ordered_samples=[{
            "GT": "0/1",
            "AD": ".,.",
            "UC": "1",
            "UV": "1",
            "LN": "10",
            "ST": ".",
            "QV": "60",
            "TY": "INS",
            "ID": source_id,
            "SC": "OctopuSV",
            "REF": "N",
            "ALT": "<INS>",
            "CO": "chr1_100-chr1_110",
        }],
    )


@pytest.mark.parametrize("unsafe", ["bad,name", "bad;name", "bad=name", "bad name"])
def test_multi_writer_rejects_unsafe_source_labels(tmp_path, unsafe):
    mapper = NameMapper(
        [str(tmp_path / "sample.svcf")],
        mode="sample",
        custom_names=[unsafe],
    )
    writer = MultiSampleWriter(mapper)

    with pytest.raises(ValueError, match="SOURCES item"):
        writer._write_event(StringIO(), _sample_mode_event())


@pytest.mark.parametrize("unsafe", ["bad,id", "bad;id", "bad=id", "bad id"])
def test_multi_writer_rejects_unsafe_source_ids(tmp_path, unsafe):
    mapper = NameMapper([str(tmp_path / "sample.svcf")], mode="sample")
    writer = MultiSampleWriter(mapper)

    with pytest.raises(ValueError, match="SOURCE_IDS item"):
        writer._write_event(StringIO(), _sample_mode_event(source_id=unsafe))


def test_validator_rejects_unsafe_versioned_source_label(tmp_path):
    path = _write_record_file(
        tmp_path / "unsafe_source.svcf",
        meta=[
            "##fileformat=VCFv4.2",
            f"##SVCFVersion={SVCF_VERSION}",
            "##OctopuSV_mode=caller",
        ],
        fmt=CALLER_FORMAT,
        blocks=[CALLER_BLOCK],
        info_suffix="SOURCES=bad=name;SOURCE_IDS=caller.1",
    )

    validator = SVCFValidator(str(path))
    validator.validate()

    assert "E_SRC_005" in {issue.code for issue in validator.errors}


def test_validator_rejects_unsafe_versioned_source_id(tmp_path):
    bad_block = CALLER_BLOCK.replace("caller.1", "bad=id")
    path = _write_record_file(
        tmp_path / "unsafe_source_id.svcf",
        meta=[
            "##fileformat=VCFv4.2",
            f"##SVCFVersion={SVCF_VERSION}",
            "##OctopuSV_mode=caller",
        ],
        fmt=CALLER_FORMAT,
        blocks=[bad_block],
        info_suffix="SOURCES=caller;SOURCE_IDS=bad=id",
    )

    validator = SVCFValidator(str(path))
    validator.validate()

    assert "E_SRC_006" in {issue.code for issue in validator.errors}


def test_legacy_validator_does_not_retroactively_apply_v11_source_atom_rule(tmp_path):
    path = _write_record_file(
        tmp_path / "legacy_unsafe_label.svcf",
        meta=["##fileformat=VCFv4.2"],
        fmt=CALLER_FORMAT,
        blocks=[CALLER_BLOCK],
        info_suffix="SOURCES=bad=name;SOURCE_IDS=caller.1",
    )

    validator = SVCFValidator(str(path))
    validator.validate()

    assert "E_SRC_005" not in {issue.code for issue in validator.errors}
