from pathlib import Path
from types import SimpleNamespace

import pytest

from octopusv.cli.merge import _preflight_merge_inputs
from octopusv.filtering.svcf_filter import FilterConfig, SVCFFilter
from octopusv.formatter.svcf_to_bedpe_converter import SVCFtoBEDPEConverter
from octopusv.formatter.svcf_to_vcf_converter import SVCFtoVCFConverter
from octopusv.utils.svcf_schema import CALLER_FORMAT
from octopusv.utils.svcf_validator import SVCFValidator


BASE_INFO = (
    "SVTYPE={svtype};END={end};SVLEN={svlen};CHR2={chr2};SUPPORT=5;"
    "SVMETHOD=OctopuSV;RTID=.;AF=.;STRAND={strand};RNAMES=."
)


def _evidence(*, svtype="DEL", record_id="r1", caller="caller", alt="<DEL>", co=None):
    if co is None:
        co = "chr1_100-chr1_200"
    ln = "." if svtype in {"TRA", "BND"} else "100"
    return (
        f"0/1:5,5:{ln}:.:60:{svtype}:{record_id}:{caller}:N:{alt}:{co}"
    )


def _write_caller_svcf(
    path: Path,
    *,
    chrom="chr1",
    pos=100,
    record_id="r1",
    alt="<DEL>",
    svtype="DEL",
    end="200",
    svlen="100",
    chr2=None,
    strand=".",
    extra_info="",
    extra_meta=None,
    contigs=None,
    versioned=True,
):
    chr2 = chrom if chr2 is None else chr2
    if contigs is None:
        contigs = [chrom] if chr2 == chrom else [chrom, chr2]

    meta = ["##fileformat=VCFv4.2"]
    if versioned:
        meta += ["##SVCFVersion=1.1", "##OctopuSV_mode=caller"]
    for contig in contigs:
        meta.append(f"##contig=<ID={contig},length=1000000>")
    if extra_meta:
        meta.extend(extra_meta)

    info = BASE_INFO.format(
        svtype=svtype,
        end=end,
        svlen=svlen,
        chr2=chr2,
        strand=strand,
    )
    if extra_info:
        info += ";" + extra_info

    co = f"{chrom}_{pos}-{chr2}_{end}" if str(end).isdigit() else "."
    block = _evidence(
        svtype=svtype,
        record_id=record_id,
        alt=alt,
        co=co,
    )
    row = (
        f"{chrom}\t{pos}\t{record_id}\tN\t{alt}\t60\tPASS\t{info}\t"
        f"{CALLER_FORMAT}\t{block}"
    )
    path.write_text(
        "\n".join(
            meta
            + [
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE",
                row,
            ]
        )
        + "\n"
    )


def _data_rows(path: Path):
    return [line for line in path.read_text().splitlines() if line and not line.startswith("#")]


def test_svcf2vcf_preserves_custom_info_flag_and_value(tmp_path):
    source = tmp_path / "input.svcf"
    output = tmp_path / "output.vcf"
    _write_caller_svcf(
        source,
        extra_meta=[
            '##INFO=<ID=PRECISE,Number=0,Type=Flag,Description="Precise breakpoint">',
            '##INFO=<ID=CUSTOM,Number=1,Type=String,Description="Caller annotation">',
        ],
        extra_info="PRECISE;CUSTOM=abc",
    )

    SVCFtoVCFConverter(None, source).convert_to_file(output)

    row_info = _data_rows(output)[0].split("\t")[7].split(";")
    assert "PRECISE" in row_info
    assert "PRECISE=True" not in row_info
    assert "CUSTOM=abc" in row_info
    text = output.read_text()
    assert "##INFO=<ID=PRECISE," in text
    assert "##INFO=<ID=CUSTOM," in text


def test_svcf2vcf_preserves_custom_alt_definition(tmp_path):
    source = tmp_path / "input.svcf"
    output = tmp_path / "output.vcf"
    custom_alt = "<DEL:ME:ALU>"
    _write_caller_svcf(
        source,
        alt=custom_alt,
        extra_meta=[
            '##ALT=<ID=DEL:ME:ALU,Description="Mobile element deletion representation">'
        ],
    )

    SVCFtoVCFConverter(None, source).convert_to_file(output)

    text = output.read_text()
    assert '##ALT=<ID=DEL:ME:ALU,Description="Mobile element deletion representation">' in text
    assert _data_rows(output)[0].split("\t")[4] == custom_alt


def test_svcf2vcf_symbolic_tra_preserves_remote_breakpoint_in_chr2_end(tmp_path):
    source = tmp_path / "tra.svcf"
    output = tmp_path / "tra.vcf"
    _write_caller_svcf(
        source,
        alt="<TRA>",
        svtype="TRA",
        end="5000",
        svlen=".",
        chr2="chr5",
        contigs=["chr1", "chr5"],
    )

    SVCFtoVCFConverter(None, source).convert_to_file(output)

    fields = _data_rows(output)[0].split("\t")
    assert fields[4] == "<TRA>"
    info_items = set(fields[7].split(";"))
    assert "CHR2=chr5" in info_items
    assert "END=5000" in info_items


def _write_filter_input(path: Path, records):
    header = (
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n"
    )
    rows = []
    for idx, (record_id, sources, sc) in enumerate(records, start=1):
        info = (
            f"SVTYPE=DEL;END={100 + idx};SVLEN=1;CHR2=chr1;SUPPORT=5;"
            "SVMETHOD=OctopuSV;RTID=.;AF=.;STRAND=.;RNAMES=."
        )
        if sources is not None:
            info += f";SOURCES={sources};SOURCE_IDS={record_id}"
        block = _evidence(record_id=record_id, caller=sc, co=f"chr1_{idx}-chr1_{100 + idx}")
        rows.append(
            f"chr1\t{idx}\t{record_id}\tN\t<DEL>\t60\tPASS\t{info}\t{CALLER_FORMAT}\t{block}"
        )
    path.write_text(header + "\n".join(rows) + "\n")


def test_filter_source_identity_is_case_sensitive_and_counts_distinct_labels(tmp_path):
    source = tmp_path / "sources.svcf"
    output = tmp_path / "filtered.svcf"
    _write_filter_input(source, [("mixed", "A,a", "caller")])

    summary = SVCFFilter(
        FilterConfig(
            input_file=str(source),
            output_file=str(output),
            min_source_count=2,
        )
    ).run()

    assert summary["output_records"] == 1
    assert _data_rows(output)[0].split("\t")[2] == "mixed"


@pytest.mark.parametrize(
    "sources,sc",
    [
        ("PBSV", "caller"),
        (None, "PBSV"),
    ],
)
def test_filter_wrong_case_source_gives_actionable_error(tmp_path, sources, sc):
    source = tmp_path / "source.svcf"
    _write_filter_input(source, [("r1", sources, sc)])

    with pytest.raises(ValueError, match=r"No exact source 'pbsv'.*Did you mean 'PBSV'"):
        SVCFFilter(
            FilterConfig(
                input_file=str(source),
                output_file=None,
                dry_run=True,
                sources={"pbsv"},
            )
        ).run()


def test_filter_exact_source_does_not_match_case_variant(tmp_path):
    source = tmp_path / "source.svcf"
    output = tmp_path / "filtered.svcf"
    _write_filter_input(source, [("upper", "A", "caller"), ("lower", "a", "caller")])

    summary = SVCFFilter(
        FilterConfig(
            input_file=str(source),
            output_file=str(output),
            sources={"A"},
        )
    ).run()

    assert summary["output_records"] == 1
    assert _data_rows(output)[0].split("\t")[2] == "upper"


def test_merge_versioned_v11_rejects_bad_del_end_via_shared_semantics(tmp_path):
    source = tmp_path / "bad_end.svcf"
    _write_caller_svcf(source, end="oops")

    with pytest.raises(ValueError, match=r"E_END_002.*DEL END must be an integer"):
        _preflight_merge_inputs(
            input_files=[source],
            labels=["caller"],
            mode="caller",
        )


def test_merge_versioned_v11_rejects_illegal_svtype_via_shared_semantics(tmp_path):
    source = tmp_path / "illegal_type.svcf"
    _write_caller_svcf(
        source,
        alt="<CPX>",
        svtype="CPX",
    )

    with pytest.raises(ValueError, match=r"E_SVTYPE_001.*Illegal SVTYPE 'CPX'"):
        _preflight_merge_inputs(
            input_files=[source],
            labels=["caller"],
            mode="caller",
        )


def test_merge_versioned_v11_rejects_symbolic_bnd_via_shared_semantics(tmp_path):
    source = tmp_path / "symbolic_bnd.svcf"
    _write_caller_svcf(
        source,
        alt="<BND>",
        svtype="BND",
        end="5000",
        svlen=".",
        chr2="chr5",
        contigs=["chr1", "chr5"],
    )

    with pytest.raises(ValueError, match=r"E_TRA_003.*BND ALT must be a valid BND bracket form"):
        _preflight_merge_inputs(
            input_files=[source],
            labels=["caller"],
            mode="caller",
        )


def test_merge_legacy_input_keeps_legacy_semantic_compatibility(tmp_path):
    source = tmp_path / "legacy_bad_end.svcf"
    _write_caller_svcf(source, end="oops", versioned=False)

    _preflight_merge_inputs(
        input_files=[source],
        labels=["caller"],
        mode="caller",
    )


def test_merge_accepts_versioned_symbolic_tra_contract(tmp_path):
    source = tmp_path / "tra.svcf"
    _write_caller_svcf(
        source,
        alt="<TRA>",
        svtype="TRA",
        end="5000",
        svlen=".",
        chr2="chr5",
        contigs=["chr1", "chr5"],
    )

    _preflight_merge_inputs(
        input_files=[source],
        labels=["caller"],
        mode="caller",
    )


def test_merge_rejects_definite_chr_prefix_alias_conflict(tmp_path):
    prefixed = tmp_path / "prefixed.svcf"
    bare = tmp_path / "bare.svcf"
    _write_caller_svcf(prefixed, chrom="chr1", contigs=["chr1"])
    _write_caller_svcf(bare, chrom="1", chr2="1", contigs=["1"])

    with pytest.raises(ValueError, match=r"contig naming.*chr1.*1|chr1.*1.*contig naming"):
        _preflight_merge_inputs(
            input_files=[prefixed, bare],
            labels=["prefixed", "bare"],
            mode="caller",
        )


def test_merge_does_not_guess_mismatch_for_disjoint_standard_contigs(tmp_path):
    chr1 = tmp_path / "chr1.svcf"
    bare2 = tmp_path / "bare2.svcf"
    _write_caller_svcf(chr1, chrom="chr1", contigs=["chr1"])
    _write_caller_svcf(bare2, chrom="2", chr2="2", contigs=["2"])

    _preflight_merge_inputs(
        input_files=[chr1, bare2],
        labels=["a", "b"],
        mode="caller",
    )


@pytest.mark.parametrize("svtype", ["DEL", "TRA"])
def test_bedpe_unknown_strand_stays_unknown(tmp_path, svtype):
    event = SimpleNamespace(
        chrom="chr1",
        pos=100,
        sv_id="event1",
        sv_type=svtype,
        quality="60",
        info={
            "CHR2": "chr5" if svtype == "TRA" else "chr1",
            "END": "5000" if svtype == "TRA" else "200",
            "SVLEN": "." if svtype == "TRA" else "100",
            "SUPPORT": "5",
            "STRAND": ".",
        },
    )

    row = SVCFtoBEDPEConverter([event]).convert().splitlines()[1].split("\t")
    assert row[8:10] == [".", "."]


def test_bedpe_known_strand_is_preserved():
    event = SimpleNamespace(
        chrom="chr1",
        pos=100,
        sv_id="event1",
        sv_type="TRA",
        quality="60",
        info={
            "CHR2": "chr5",
            "END": "5000",
            "SVLEN": ".",
            "SUPPORT": "5",
            "STRAND": "+-",
        },
    )

    row = SVCFtoBEDPEConverter([event]).convert().splitlines()[1].split("\t")
    assert row[8:10] == ["+", "-"]


@pytest.mark.parametrize(
    "chrom,chr2",
    [
        ("HLA-A*01:01:01:01", "HLA-A*01:01:01:01"),
        ("chr1", "HLA-A*01:01:01:01"),
    ],
)
def test_validator_rejects_colon_contigs_that_cannot_round_trip_svcf_evidence(
    tmp_path, chrom, chr2
):
    source = tmp_path / "colon_contig.svcf"
    _write_caller_svcf(
        source,
        chrom=chrom,
        chr2=chr2,
        contigs=[chrom] if chrom == chr2 else [chrom, chr2],
    )

    validator = SVCFValidator(str(source))
    validator.validate()

    assert any(issue.code == "E_CONTIG_001" for issue in validator.errors), (
        validator.to_summary(max_issues=None)
    )


def test_merge_preflight_reuses_colon_contig_contract_check(tmp_path):
    source = tmp_path / "colon_contig.svcf"
    _write_caller_svcf(
        source,
        chrom="HLA-A*01:01:01:01",
        chr2="HLA-A*01:01:01:01",
        contigs=["HLA-A*01:01:01:01"],
    )

    with pytest.raises(ValueError, match=r"E_CONTIG_001.*contig.*:"):
        _preflight_merge_inputs(
            input_files=[source],
            labels=["caller"],
            mode="caller",
        )


def test_merge_preflight_rejects_local_colon_contig_with_normal_remote_chr2(tmp_path):
    """Guard the local-CHROM branch of E_CONTIG_001 independently of CHR2."""
    source = tmp_path / "local_colon_remote_normal.svcf"
    _write_caller_svcf(
        source,
        chrom="HLA-A*01:01:01:01",
        chr2="chr5",
        alt="<TRA>",
        svtype="TRA",
        end="5000",
        svlen=".",
        strand=".",
        contigs=["HLA-A*01:01:01:01", "chr5"],
    )

    with pytest.raises(ValueError, match=r"E_CONTIG_001.*CHROM=.*:"):
        _preflight_merge_inputs(
            input_files=[source],
            labels=["caller"],
            mode="caller",
        )
