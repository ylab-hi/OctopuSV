from __future__ import annotations

import builtins
from pathlib import Path

import pytest

from octopusv.utils.svcf_validator import SVCFValidator


FORMAT = "GT:AD:LN:ST:QV:TY:ID:SC:REF:ALT:CO"
BASE_INFO = (
    "SVTYPE=INS;END=110;SVLEN=10;CHR2=chr1;SUPPORT=5;"
    "SVMETHOD=OctopuSV;RTID=.;AF=.;STRAND=.;RNAMES=."
)


def _evidence(source_id: str, sc: str = "caller", co: str = "chr1_100-chr1_110") -> str:
    return (
        "0/1:5,7:10:.:60:INS:"
        f"{source_id}:{sc}:N:<INS>:{co}"
    )


def _placeholder() -> str:
    return "0/0:.:.:.:.:.:.:.:.:.:."


def _write_svcf(
    path: Path,
    *,
    info_suffix: str = "",
    evidence: list[str] | None = None,
    multi: bool = False,
    sample_names: list[str] | None = None,
    info_override: str | None = None,
    alt: str = "<INS>",
) -> Path:
    if evidence is None:
        evidence = [_evidence("source.1")]

    header = ["##fileformat=VCFv4.2"]
    if multi:
        header.append("##OctopuSV_mode=multi")

    if sample_names is None:
        sample_names = ["SAMPLE"]

    header.append(
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"
        + "\t".join(sample_names)
    )

    if info_override is not None:
        info = info_override
    else:
        info = BASE_INFO + (";" + info_suffix if info_suffix else "")

    record = (
        "chr1\t100\tmerged.1\tN\t"
        f"{alt}\t60\tPASS\t{info}\t{FORMAT}\t"
        + "\t".join(evidence)
    )

    path.write_text("\n".join(header + [record]) + "\n", encoding="utf-8")
    return path


def _validate(path: Path, *, strict_co: bool = False) -> SVCFValidator:
    validator = SVCFValidator(str(path), strict_co=strict_co)
    validator.validate()
    return validator


def _codes(validator: SVCFValidator) -> set[str]:
    return {issue.code for issue in validator.issues}


def test_caller_merge_allows_duplicate_sources_when_ids_and_evidence_align(tmp_path):
    path = _write_svcf(
        tmp_path / "duplicate_sources.svcf",
        info_suffix=(
            "SOURCES=sniffles,sniffles,pbsv;"
            "SOURCE_IDS=sniffles.1,sniffles.2,pbsv.1"
        ),
        evidence=[
            _evidence("sniffles.1", sc="Sniffles2-2.2"),
            _evidence("sniffles.2", sc="Sniffles2-2.2"),
            _evidence("pbsv.1", sc="pbsv-2.9"),
        ],
    )

    validator = _validate(path)

    assert validator.mode == "caller_merge"
    assert validator.status() == "PASS"


def test_source_ids_remain_optional_for_legacy_caller_merge(tmp_path):
    path = _write_svcf(
        tmp_path / "legacy.svcf",
        info_suffix="SOURCES=sourceA,sourceB",
        evidence=[_evidence("idA"), _evidence("idB")],
    )

    validator = _validate(path)
    assert validator.errors == []


def test_source_ids_count_must_match_sources_when_present(tmp_path):
    path = _write_svcf(
        tmp_path / "bad_count.svcf",
        info_suffix="SOURCES=sourceA,sourceB;SOURCE_IDS=idA",
        evidence=[_evidence("idA"), _evidence("idB")],
    )

    assert "E_SRC_003" in _codes(_validate(path))


def test_source_ids_match_evidence_ids_by_position(tmp_path):
    path = _write_svcf(
        tmp_path / "bad_order.svcf",
        info_suffix="SOURCES=sourceA,sourceB;SOURCE_IDS=idA,idB",
        evidence=[_evidence("idB"), _evidence("idA")],
    )

    validator = _validate(path)
    assert len([x for x in validator.issues if x.code == "E_SRC_004"]) == 2


def test_source_ids_dot_placeholder_keeps_its_position(tmp_path):
    path = _write_svcf(
        tmp_path / "missing_middle_id.svcf",
        info_suffix="SOURCES=a,b,c;SOURCE_IDS=idA,.,idC",
        evidence=[_evidence("idA"), _evidence("."), _evidence("idC")],
    )

    validator = _validate(path)
    assert validator.errors == []


def test_colon_containing_source_id_uses_shared_parser(tmp_path):
    manta_id = "MantaDEL:469174:0:1:0:0:0"
    path = _write_svcf(
        tmp_path / "manta.svcf",
        info_suffix=f"SOURCES=manta;SOURCE_IDS={manta_id}",
        evidence=[_evidence(manta_id, sc="Manta_v1.6.0")],
    )

    assert _validate(path).errors == []


def test_validator_does_not_guess_identity_from_sc(tmp_path):
    path = _write_svcf(
        tmp_path / "no_guess.svcf",
        info_suffix="SOURCES=user_label;SOURCE_IDS=record.1",
        evidence=[_evidence("record.1", sc="different_name")],
    )

    assert _validate(path).errors == []


def test_sample_multi_only_checks_source_id_list_count(tmp_path):
    path = _write_svcf(
        tmp_path / "multi.svcf",
        info_suffix="SOURCES=sampleA,sampleC;SOURCE_IDS=idA,idC",
        evidence=[_evidence("idA"), _placeholder(), _evidence("idC")],
        multi=True,
        sample_names=["sampleA", "sampleB", "sampleC"],
    )

    validator = _validate(path)
    assert validator.mode == "sample_multi"
    assert validator.errors == []


def test_duplicate_info_keys_are_rejected(tmp_path):
    path = _write_svcf(
        tmp_path / "duplicate_info.svcf",
        info_override=BASE_INFO + ";SOURCES=a;SOURCE_IDS=id1;SOURCE_IDS=id1",
        evidence=[_evidence("id1")],
    )

    assert "E_INFO_002" in _codes(_validate(path))


def test_truncated_fixed_evidence_block_is_rejected(tmp_path):
    truncated = "0/1:5,7:10:.:60:INS:id1:caller:N"
    path = _write_svcf(
        tmp_path / "truncated.svcf",
        info_suffix="SOURCES=a;SOURCE_IDS=id1",
        evidence=[truncated],
    )

    assert "E_FMT_002" in _codes(_validate(path))


def test_co_with_hyphenated_contigs_is_valid(tmp_path):
    path = _write_svcf(
        tmp_path / "hyphen_co.svcf",
        evidence=[_evidence("id1", co="chr1-alt_100-chr2-alt_110")],
    )

    validator = _validate(path, strict_co=True)
    assert "E_CO_001" not in _codes(validator)


def test_ambiguous_co_warns_or_errors_in_strict_mode(tmp_path):
    ambiguous = "chrY_KI270740v1_random_123-NC_007605-alt_456"
    path = _write_svcf(
        tmp_path / "ambiguous_co.svcf",
        evidence=[_evidence("id1", co=ambiguous)],
    )

    assert "W_CO_001" in _codes(_validate(path))
    assert "E_CO_001" in _codes(_validate(path, strict_co=True))


def test_mixed_no_marker_mode_is_rejected_streaming(tmp_path):
    path = tmp_path / "mixed.svcf"
    record1 = (
        f"chr1\t100\tr1\tN\t<INS>\t60\tPASS\t{BASE_INFO}\t{FORMAT}\t"
        + _evidence("id1")
    )
    record2 = (
        f"chr1\t200\tr2\tN\t<INS>\t60\tPASS\t{BASE_INFO};"
        f"SOURCES=a;SOURCE_IDS=id2\t{FORMAT}\t"
        + _evidence("id2", co="chr1_200-chr1_210")
    )
    path.write_text(
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n"
        + record1
        + "\n"
        + record2
        + "\n",
        encoding="utf-8",
    )

    assert "E_MODE_002" in _codes(_validate(path))


def test_validator_streams_and_never_calls_readlines(tmp_path, monkeypatch):
    path = _write_svcf(tmp_path / "streaming.svcf")
    real_open = builtins.open

    class NoReadlinesWrapper:
        def __init__(self, handle):
            self._handle = handle

        def __iter__(self):
            return iter(self._handle)

        def __enter__(self):
            self._handle.__enter__()
            return self

        def __exit__(self, *args):
            return self._handle.__exit__(*args)

        def readlines(self, *args, **kwargs):
            raise AssertionError("validator must stream instead of calling readlines()")

    def wrapped_open(file, *args, **kwargs):
        handle = real_open(file, *args, **kwargs)
        if str(file) == str(path):
            return NoReadlinesWrapper(handle)
        return handle

    monkeypatch.setattr(builtins, "open", wrapped_open)
    assert _validate(path).status() in {"PASS", "PASS_WITH_WARNINGS"}


def test_tra_requires_matching_brackets_and_dot_svlen(tmp_path):
    info = (
        "SVTYPE=TRA;END=200;SVLEN=10;CHR2=chr2;SUPPORT=5;"
        "SVMETHOD=OctopuSV;RTID=.;AF=.;STRAND=.;RNAMES=."
    )
    path = _write_svcf(
        tmp_path / "bad_tra.svcf",
        info_override=info,
        alt="N[chr2:200]",
        evidence=[_evidence("id1")],
    )

    codes = _codes(_validate(path))
    assert "E_TRA_003" in codes
    assert "E_TRA_006" in codes


def test_summary_gives_actionable_migration_note_for_source_or_format_errors(tmp_path):
    path = _write_svcf(
        tmp_path / "old_affected.svcf",
        info_suffix="SOURCES=a;SOURCE_IDS=idA",
        evidence=[_evidence("different_id")],
    )

    validator = _validate(path)
    summary = validator.to_summary()

    assert "E_SRC_004" in _codes(validator)
    assert "Migration note:" in summary
    assert "re-run the original merge with OctopuSV 0.5.0" in summary
    assert "manually editing SOURCE_IDS" in summary
