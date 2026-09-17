from __future__ import annotations

from pathlib import Path

import pytest
from typer.testing import CliRunner

pytest.importorskip("natsort")

from octopusv.cli.cli import app
from octopusv.normalization.contig_normalizer import SVCFContigNormalizer


runner = CliRunner()
CALLER_FORMAT = "GT:AD:LN:ST:QV:TY:ID:SC:REF:ALT:CO"


def _invoke(args: list[str]):
    result = runner.invoke(app, args)
    assert result.exit_code == 0, (
        f"Command failed: octopusv {' '.join(args)}\n"
        f"output:\n{result.output}\n"
        f"exception: {result.exception!r}"
    )
    return result


def _meta_lines(path: Path) -> list[str]:
    return [
        line
        for line in path.read_text(encoding="utf-8").splitlines()
        if line.startswith("##")
    ]


def _write_raw_del(path: Path, *, reference: str, assembly: str) -> None:
    path.write_text(
        "\n".join(
            [
                "##fileformat=VCFv4.2",
                f"##reference={reference}",
                f"##assembly={assembly}",
                "##source=testcaller",
                "##contig=<ID=chr1,length=1000000>",
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1",
                "chr1\t100\traw1\tN\t<DEL>\t60\tPASS\tSVTYPE=DEL;END=600;SVLEN=-500;RE=10\tGT\t0/1",
            ]
        )
        + "\n",
        encoding="utf-8",
    )


def _caller_block(record_id: str, source: str, *, pos: int, end: int) -> str:
    return ":".join(
        [
            "0/1",
            ".,.",
            str(abs(end - pos)),
            ".",
            "60",
            "DEL",
            record_id,
            source,
            "N",
            "<DEL>",
            f"chr1_{pos}-chr1_{end}",
        ]
    )


def _write_caller_svcf(
    path: Path,
    *,
    record_id: str,
    pos: int,
    source: str,
    reference: str | None = "GRCh38",
    assembly: str | None = "GRCh38.p14",
) -> None:
    end = pos + 500
    meta = [
        "##fileformat=VCFv4.2",
        "##SVCFVersion=1.1",
        "##OctopuSV_mode=caller",
    ]
    if reference is not None:
        meta.append(f"##reference={reference}")
    if assembly is not None:
        meta.append(f"##assembly={assembly}")
    meta += [
        "##contig=<ID=chr1,length=1000000>",
        '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="SV type">',
        '##INFO=<ID=END,Number=1,Type=Integer,Description="End">',
        '##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Length">',
        '##INFO=<ID=CHR2,Number=1,Type=String,Description="Mate contig">',
        '##INFO=<ID=SUPPORT,Number=1,Type=Integer,Description="Support">',
        '##INFO=<ID=SVMETHOD,Number=1,Type=String,Description="Method">',
        '##INFO=<ID=RTID,Number=1,Type=String,Description="Related ID">',
        '##INFO=<ID=AF,Number=1,Type=Float,Description="AF">',
        '##INFO=<ID=STRAND,Number=1,Type=String,Description="Strand">',
        '##INFO=<ID=RNAMES,Number=.,Type=String,Description="Reads">',
        '##INFO=<ID=SOURCES,Number=.,Type=String,Description="Sources">',
        '##INFO=<ID=SOURCE_IDS,Number=.,Type=String,Description="Source IDs">',
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
        '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depth">',
        '##FORMAT=<ID=LN,Number=1,Type=Integer,Description="Length">',
        '##FORMAT=<ID=ST,Number=1,Type=String,Description="Strand">',
        '##FORMAT=<ID=QV,Number=1,Type=Float,Description="Quality">',
        '##FORMAT=<ID=TY,Number=1,Type=String,Description="Type">',
        '##FORMAT=<ID=ID,Number=1,Type=String,Description="ID">',
        '##FORMAT=<ID=SC,Number=1,Type=String,Description="Source caller">',
        '##FORMAT=<ID=REF,Number=1,Type=String,Description="REF">',
        '##FORMAT=<ID=ALT,Number=1,Type=String,Description="ALT">',
        '##FORMAT=<ID=CO,Number=1,Type=String,Description="Coordinate">',
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE",
    ]
    info = (
        f"SVTYPE=DEL;END={end};SVLEN=500;CHR2=chr1;SUPPORT=10;"
        "SVMETHOD=OctopuSV;RTID=.;AF=.;STRAND=.;RNAMES=.;"
        f"SOURCES={source};SOURCE_IDS={record_id}"
    )
    block = _caller_block(record_id, source, pos=pos, end=end)
    path.write_text(
        "\n".join(
            meta
            + [
                f"chr1\t{pos}\t{record_id}\tN\t<DEL>\t60\tPASS\t{info}\t{CALLER_FORMAT}\t{block}"
            ]
        )
        + "\n",
        encoding="utf-8",
    )


def test_reference_and_assembly_survive_correct_and_svcf2vcf(tmp_path):
    raw = tmp_path / "input.vcf"
    corrected = tmp_path / "corrected.svcf"
    exported = tmp_path / "exported.vcf"
    _write_raw_del(raw, reference="GRCh38", assembly="GRCh38.p14")

    _invoke(["correct", "-i", str(raw), "-o", str(corrected)])
    corrected_meta = _meta_lines(corrected)
    assert corrected_meta.count("##reference=GRCh38") == 1
    assert corrected_meta.count("##assembly=GRCh38.p14") == 1

    _invoke(["validate-svcf", str(corrected)])
    _invoke(["svcf2vcf", "-i", str(corrected), "-o", str(exported)])
    exported_meta = _meta_lines(exported)
    assert exported_meta.count("##reference=GRCh38") == 1
    assert exported_meta.count("##assembly=GRCh38.p14") == 1


@pytest.mark.parametrize("mode", ["caller", "sample"])
def test_merge_and_export_preserve_consistent_global_metadata(tmp_path, mode):
    a = tmp_path / "a.svcf"
    b = tmp_path / "b.svcf"
    merged = tmp_path / f"merged.{mode}.svcf"
    exported = tmp_path / f"merged.{mode}.vcf"
    _write_caller_svcf(a, record_id="a1", pos=100, source="callerA")
    # Missing metadata in one input is not an explicit conflict.
    _write_caller_svcf(
        b,
        record_id="b1",
        pos=105,
        source="callerB",
        reference=None,
        assembly=None,
    )

    args = [
        "merge",
        "-i", str(a),
        "-i", str(b),
        "-o", str(merged),
        "--mode", mode,
        "--union",
    ]
    if mode == "caller":
        args += ["--caller-names", "callerA,callerB"]
    else:
        args += ["--sample-names", "sampleA,sampleB"]
    _invoke(args)

    merged_meta = _meta_lines(merged)
    assert merged_meta.count("##reference=GRCh38") == 1
    assert merged_meta.count("##assembly=GRCh38.p14") == 1

    _invoke(["validate-svcf", str(merged)])
    _invoke(["svcf2vcf", "-i", str(merged), "-o", str(exported)])
    exported_meta = _meta_lines(exported)
    assert exported_meta.count("##reference=GRCh38") == 1
    assert exported_meta.count("##assembly=GRCh38.p14") == 1


@pytest.mark.parametrize("mode", ["caller", "sample"])
@pytest.mark.parametrize(
    ("field", "left", "right"),
    [
        ("reference", "GRCh37", "GRCh38"),
        ("assembly", "GRCh38.p13", "GRCh38.p14"),
    ],
)
def test_merge_rejects_explicit_global_metadata_conflicts(
    tmp_path, mode, field, left, right
):
    a = tmp_path / "a.svcf"
    b = tmp_path / "b.svcf"
    out = tmp_path / "out.svcf"

    kwargs_a = {"reference": "GRCh38", "assembly": "GRCh38.p14"}
    kwargs_b = dict(kwargs_a)
    kwargs_a[field] = left
    kwargs_b[field] = right
    _write_caller_svcf(a, record_id="a1", pos=100, source="callerA", **kwargs_a)
    _write_caller_svcf(b, record_id="b1", pos=105, source="callerB", **kwargs_b)

    args = [
        "merge",
        "-i", str(a),
        "-i", str(b),
        "-o", str(out),
        "--mode", mode,
        "--union",
    ]
    if mode == "caller":
        args += ["--caller-names", "callerA,callerB"]
    else:
        args += ["--sample-names", "sampleA,sampleB"]

    result = runner.invoke(app, args)
    assert result.exit_code != 0
    assert f"Conflicting ##{field} metadata" in result.output or (
        result.exception is not None
        and f"Conflicting ##{field} metadata" in str(result.exception)
    )
    assert not out.exists()


def _minimal_normalize_svcf(
    *,
    contigs: list[str],
    chrom: str,
    chr2: str,
    alt: str = "<DEL>",
    co: str | None = None,
) -> str:
    if co is None:
        co = f"{chrom}_100-{chr2}_600"
    meta = [
        "##fileformat=VCFv4.2",
        "##SVCFVersion=1.1",
        "##OctopuSV_mode=caller",
        *[f"##contig=<ID={contig},length=1000000>" for contig in contigs],
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE",
    ]
    info = (
        f"SVTYPE=DEL;END=600;SVLEN=500;CHR2={chr2};SUPPORT=10;"
        "SVMETHOD=OctopuSV;RTID=.;AF=.;STRAND=.;RNAMES=."
    )
    block = f"0/1:.,.:500:.:60:DEL:r1:test:N:{alt}:{co}"
    return (
        "\n".join(meta)
        + "\n"
        + f"{chrom}\t100\tr1\tN\t{alt}\t60\tPASS\t{info}\t{CALLER_FORMAT}\t{block}\n"
    )


def test_normalize_contigs_rejects_header_many_to_one_collision_atomically(tmp_path):
    source = tmp_path / "mixed.svcf"
    output = tmp_path / "out.svcf"
    source.write_text(
        _minimal_normalize_svcf(
            contigs=["1", "chr1"],
            chrom="1",
            chr2="1",
            co="1_100-1_600",
        ),
        encoding="utf-8",
    )
    output.write_text("ORIGINAL\n", encoding="utf-8")

    with pytest.raises(ValueError, match=r"would both normalize to 'chr1'"):
        SVCFContigNormalizer(source, output, "chr").run()

    assert output.read_text(encoding="utf-8") == "ORIGINAL\n"
    assert not list(tmp_path.glob(f".{output.name}.*.tmp"))


@pytest.mark.parametrize(
    ("contigs", "chrom", "chr2", "alt", "co"),
    [
        # Header uses chr1, record CHROM uses 1.
        (["chr1"], "1", "1", "<DEL>", "1_100-1_600"),
        # CHROM/header use 1, INFO/CHR2 switches to chr1.
        (["1"], "1", "chr1", "<DEL>", "1_100-1_600"),
        # CHROM/header use 1, BND ALT mate switches to chr1.
        (["1"], "1", "1", "N]chr1:600]", "1_100-1_600"),
        # CHROM/header use 1, FORMAT/CO switches one endpoint to chr1.
        (["1"], "1", "1", "<DEL>", "1_100-chr1_600"),
    ],
    ids=["record_chrom", "info_chr2", "bnd_alt", "format_co"],
)
def test_normalize_contigs_rejects_cross_field_namespace_collapse(
    tmp_path, contigs, chrom, chr2, alt, co
):
    source = tmp_path / "mixed.svcf"
    output = tmp_path / "out.svcf"
    source.write_text(
        _minimal_normalize_svcf(
            contigs=contigs,
            chrom=chrom,
            chr2=chr2,
            alt=alt,
            co=co,
        ),
        encoding="utf-8",
    )

    with pytest.raises(ValueError, match=r"would both normalize to 'chr1'"):
        SVCFContigNormalizer(source, output, "chr").run()

    assert not output.exists()
