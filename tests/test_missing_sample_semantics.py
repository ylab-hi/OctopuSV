from pathlib import Path
from types import SimpleNamespace

import typer
from typer.testing import CliRunner

from octopusv.cli.svcf2vcf import svcf2vcf
from octopusv.formatter.svcf_to_vcf_converter import SVCFtoVCFConverter
from octopusv.merger.multi_sample_writer import MultiSampleWriter
from octopusv.merger.name_mapper import NameMapper
from octopusv.stater.stat_analyzers import GenotypeAnalyzer
from octopusv.stater.stat_reader import read_records
from octopusv.stater.sv_stater import SVStater
from octopusv.utils.sample_mode_semantics import (
    downstream_sample_gt,
    is_unobserved_sample,
    normalize_unobserved_sample_gt,
)


FORMAT_V11 = "GT:AD:UC:UV:LN:ST:QV:TY:ID:SC:REF:ALT:CO"
FORMAT_LEGACY = "GT:AD:LN:ST:QV:TY:ID:SC:REF:ALT:CO"


def _svcf2vcf_app():
    app = typer.Typer()
    app.command(name="svcf2vcf")(svcf2vcf)
    return app


def _sample_call(gt="0/1", uc="2", uv="2", ln="500", record_id="sampleA.event"):
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
        "CO": "chr1_100-chr1_600",
    }


def _write_real_sample_mode_with_missing_sample(path: Path):
    input_files = [
        path.parent / "sampleA.caller.svcf",
        path.parent / "sampleB.caller.svcf",
    ]
    mapper = NameMapper(
        [str(p) for p in input_files],
        mode="sample",
        custom_names=["sampleA", "sampleB"],
    )
    writer = MultiSampleWriter(mapper)

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
        ordered_samples=[_sample_call(), None],
    )
    writer.write_results(path, [event], {"chr1": 1_000_000})


def _write_sample_svcf(path: Path, fmt: str, block: str):
    path.write_text(
        "##fileformat=VCFv4.2\n"
        "##OctopuSV_mode=multi\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n"
        "chr1\t100\te1\tN\t<DEL>\t60\tPASS\t"
        "SVTYPE=DEL;END=600;SVLEN=500\t"
        f"{fmt}\t{block}\n"
    )


def _vcf_data_fields(path: Path):
    return next(
        line for line in path.read_text().splitlines()
        if line and not line.startswith("#")
    ).split("\t")


def test_unobserved_placeholder_detection_is_stricter_than_uv_zero():
    placeholder = {
        "GT": "0/0", "UC": "0", "UV": "0", "ID": ".", "SC": "."
    }
    assert is_unobserved_sample(placeholder)
    assert downstream_sample_gt(placeholder) == "./."
    assert downstream_sample_gt(
        placeholder, unobserved_sample_gt="ref"
    ) == "0/0"

    # Evidence exists, but no valid vote survives.  This is genuinely
    # unresolved and must NEVER become hom-ref under the ref export policy.
    unresolved_with_evidence = {
        "GT": "./.",
        "UC": "0",
        "UV": "0",
        "ID": "sampleA.event",
        "SC": "OctopuSV",
    }
    assert not is_unobserved_sample(unresolved_with_evidence)
    assert downstream_sample_gt(
        unresolved_with_evidence, unobserved_sample_gt="ref"
    ) == "./."

    # Evidence-backed absence remains a real hom-ref call.
    true_absence = {
        "GT": "0/0", "UC": "0", "UV": "2", "ID": "e1", "SC": "OctopuSV"
    }
    assert not is_unobserved_sample(true_absence)
    assert downstream_sample_gt(true_absence) == "0/0"

    # Legacy sample blocks have no UC/UV; do not invent new semantics.
    assert not is_unobserved_sample({"GT": "0/0"})
    assert downstream_sample_gt(
        {"GT": "0/0"}, unobserved_sample_gt="ref"
    ) == "0/0"


def test_unobserved_sample_policy_validation():
    assert normalize_unobserved_sample_gt("missing") == "missing"
    assert normalize_unobserved_sample_gt("REF") == "ref"

    try:
        normalize_unobserved_sample_gt("maybe")
    except ValueError as exc:
        assert "Expected one of: missing, ref" in str(exc)
    else:
        raise AssertionError("invalid policy should raise ValueError")


def test_real_placeholder_default_exports_missing_and_stat_tracks_no_evidence(tmp_path):
    svcf = tmp_path / "population.svcf"
    vcf = tmp_path / "population.vcf"
    _write_real_sample_mode_with_missing_sample(svcf)

    # Lock the SVCF fact-layer contract: writer still emits the fixed-width
    # 0/0 placeholder.  Only downstream interpretation is configurable.
    svcf_fields = _vcf_data_fields(svcf)
    assert svcf_fields[8] == FORMAT_V11
    assert svcf_fields[9].startswith("0/1:.,.:2:2:500:")
    assert svcf_fields[10].startswith("0/0:.,.:0:0:.:")

    SVCFtoVCFConverter(None, svcf).convert_to_file(vcf)
    vcf_fields = _vcf_data_fields(vcf)
    assert vcf_fields[8] == "GT:AD:DP:UC:UV:LN"
    assert vcf_fields[9] == "0/1:.,.:.:2:2:500"
    assert vcf_fields[10] == "./.:.,.:.:0:0:."

    records, sample_names, mode = read_records(svcf)
    stats = GenotypeAnalyzer(records, sample_names, mode=mode).analyze()
    assert mode == "sample"
    assert stats["per_sample"] == {
        "sampleA": {"0/1": 1},
        "sampleB": {"./.": 1},
    }
    assert stats["overall"] == {"0/1": 1, "./.": 1}
    assert stats["no_evidence"] == {
        "overall": 1,
        "per_sample": {"sampleA": 0, "sampleB": 1},
    }


def test_ref_policy_only_changes_true_placeholder(tmp_path):
    svcf = tmp_path / "population.svcf"
    vcf = tmp_path / "population.ref.vcf"
    _write_real_sample_mode_with_missing_sample(svcf)

    SVCFtoVCFConverter(
        None,
        svcf,
        unobserved_sample_gt="ref",
    ).convert_to_file(vcf)

    fields = _vcf_data_fields(vcf)
    assert fields[9] == "0/1:.,.:.:2:2:500"
    assert fields[10] == "0/0:.,.:.:0:0:."


def test_cli_exposes_explicit_missing_vs_ref_policy(tmp_path):
    runner = CliRunner()
    app = _svcf2vcf_app()
    svcf = tmp_path / "population.svcf"
    missing_vcf = tmp_path / "missing.vcf"
    ref_vcf = tmp_path / "ref.vcf"
    _write_real_sample_mode_with_missing_sample(svcf)

    default_result = runner.invoke(
        app,
        ["-i", str(svcf), "-o", str(missing_vcf)],
    )
    assert default_result.exit_code == 0, default_result.output
    assert _vcf_data_fields(missing_vcf)[10].startswith("./.:")

    ref_result = runner.invoke(
        app,
        [
            "-i",
            str(svcf),
            "-o",
            str(ref_vcf),
            "--unobserved-sample-gt",
            "ref",
        ],
    )
    assert ref_result.exit_code == 0, ref_result.output
    assert _vcf_data_fields(ref_vcf)[10].startswith("0/0:")

    invalid_result = runner.invoke(
        app,
        [
            "-i",
            str(svcf),
            "-o",
            str(tmp_path / "bad.vcf"),
            "--unobserved-sample-gt",
            "maybe",
        ],
    )
    assert invalid_result.exit_code == 1
    assert "Expected one of: missing, ref" in invalid_result.output


def test_ref_policy_does_not_rewrite_unresolved_call_with_evidence(tmp_path):
    svcf = tmp_path / "unresolved.svcf"
    vcf = tmp_path / "unresolved.vcf"
    block = (
        "./.:.,.:0:0:500:.:60:DEL:e1:OctopuSV:"
        "N:<DEL>:chr1_100-chr1_600"
    )
    _write_sample_svcf(svcf, FORMAT_V11, block)

    SVCFtoVCFConverter(
        None,
        svcf,
        unobserved_sample_gt="ref",
    ).convert_to_file(vcf)
    fields = _vcf_data_fields(vcf)
    assert fields[9] == "./.:.,.:.:0:0:500"



def test_svstater_surfaces_no_evidence_in_structured_and_text_outputs(tmp_path):
    svcf = tmp_path / "population.svcf"
    _write_real_sample_mode_with_missing_sample(svcf)

    stater = SVStater(str(svcf))
    stater.analyze()

    payload = stater.to_json_dict()
    assert payload["genotypes"]["no_evidence"] == {
        "overall": 1,
        "per_sample": {"sampleA": 0, "sampleB": 1},
    }

    text = stater.to_text()
    assert "No-evidence sample/event cells" in text
    assert "sampleB" in text
    assert "Overall" in text

def test_evidence_backed_hom_ref_with_positive_uv_remains_hom_ref(tmp_path):
    svcf = tmp_path / "true_absence.svcf"
    vcf = tmp_path / "true_absence.vcf"
    block = (
        "0/0:.,.:0:2:500:.:60:DEL:e1:OctopuSV:"
        "N:<DEL>:chr1_100-chr1_600"
    )
    _write_sample_svcf(svcf, FORMAT_V11, block)

    SVCFtoVCFConverter(None, svcf).convert_to_file(vcf)
    fields = _vcf_data_fields(vcf)
    assert fields[9] == "0/0:.,.:.:0:2:500"

    records, sample_names, mode = read_records(svcf)
    stats = GenotypeAnalyzer(records, sample_names, mode=mode).analyze()
    assert stats["per_sample"] == {"S1": {"0/0": 1}}
    assert stats["no_evidence"] == {
        "overall": 0,
        "per_sample": {"S1": 0},
    }


def test_legacy_sample_without_uv_is_not_reinterpreted(tmp_path):
    svcf = tmp_path / "legacy.svcf"
    vcf = tmp_path / "legacy.vcf"
    block = "0/0:12,0:500:.:60:DEL:e1:callerA:N:<DEL>:chr1_100-chr1_600"
    _write_sample_svcf(svcf, FORMAT_LEGACY, block)

    SVCFtoVCFConverter(
        None,
        svcf,
        unobserved_sample_gt="ref",
    ).convert_to_file(vcf)
    fields = _vcf_data_fields(vcf)
    assert fields[9] == "0/0:12,0:12:.:.:500"

    records, sample_names, mode = read_records(svcf)
    stats = GenotypeAnalyzer(records, sample_names, mode=mode).analyze()
    assert stats["per_sample"] == {"S1": {"0/0": 1}}
    assert stats["no_evidence"] == {
        "overall": 0,
        "per_sample": {"S1": 0},
    }
