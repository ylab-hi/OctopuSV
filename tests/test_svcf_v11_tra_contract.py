from __future__ import annotations

from pathlib import Path

import pytest

from octopusv.utils.svcf_schema import CALLER_FORMAT, SVCF_VERSION
from octopusv.utils.svcf_validator import SVCFValidator


def _evidence(
    *,
    svtype: str,
    alt: str,
    co: str = "chr1_100-chr5_5000",
) -> str:
    return ":".join(
        [
            "0/1",    # GT
            ".,.",    # AD
            ".",      # LN
            ".",      # ST
            "60",     # QV
            svtype,   # TY
            "source.1",  # ID
            "caller", # SC
            "N",      # REF
            alt,       # ALT
            co,        # CO
        ]
    )


def _write_v11_event(
    path: Path,
    *,
    svtype: str,
    alt: str,
    chr2: str = "chr5",
    end: str = "5000",
    svlen: str = ".",
) -> Path:
    info = ";".join(
        [
            f"SVTYPE={svtype}",
            f"END={end}",
            f"SVLEN={svlen}",
            f"CHR2={chr2}",
            "SUPPORT=5",
            "SVMETHOD=OctopuSV",
            "RTID=.",
            "AF=.",
            "STRAND=.",
            "RNAMES=.",
        ]
    )
    co = f"chr1_100-{chr2}_{end}" if chr2 != "." and end.isdigit() else "."
    path.write_text(
        "\n".join(
            [
                "##fileformat=VCFv4.2",
                f"##SVCFVersion={SVCF_VERSION}",
                "##OctopuSV_mode=caller",
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE",
                (
                    "chr1\t100\tevent.1\tN\t"
                    f"{alt}\t60\tPASS\t{info}\t{CALLER_FORMAT}\t"
                    f"{_evidence(svtype=svtype, alt=alt, co=co)}"
                ),
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    return path


def _validate(path: Path) -> SVCFValidator:
    validator = SVCFValidator(str(path))
    validator.validate()
    return validator


def _codes(validator: SVCFValidator) -> set[str]:
    return {issue.code for issue in validator.issues}


def test_v11_symbolic_tra_with_two_known_breakpoints_is_valid(tmp_path):
    path = _write_v11_event(
        tmp_path / "symbolic_tra.svcf",
        svtype="TRA",
        alt="<TRA>",
    )

    validator = _validate(path)

    assert validator.status() == "PASS"
    assert "E_TRA_003" not in _codes(validator)


def test_v11_bracket_tra_with_matching_chr2_end_is_valid(tmp_path):
    path = _write_v11_event(
        tmp_path / "bracket_tra.svcf",
        svtype="TRA",
        alt="N]chr5:5000]",
    )

    assert _validate(path).status() == "PASS"


def test_v11_bracket_tra_rejects_conflicting_mate_chrom(tmp_path):
    path = _write_v11_event(
        tmp_path / "tra_bad_chr2.svcf",
        svtype="TRA",
        alt="N]chr6:5000]",
        chr2="chr5",
    )

    assert "E_TRA_004" in _codes(_validate(path))


def test_v11_bracket_tra_rejects_conflicting_mate_position(tmp_path):
    path = _write_v11_event(
        tmp_path / "tra_bad_end.svcf",
        svtype="TRA",
        alt="N]chr5:5001]",
        end="5000",
    )

    assert "E_TRA_005" in _codes(_validate(path))


@pytest.mark.parametrize(
    ("chr2", "end", "expected_code"),
    [
        (".", "5000", "E_TRA_001"),
        ("chr5", ".", "E_TRA_002"),
        ("chr5", "oops", "E_TRA_002"),
    ],
)
def test_v11_symbolic_tra_requires_two_known_breakpoints(
    tmp_path,
    chr2,
    end,
    expected_code,
):
    path = _write_v11_event(
        tmp_path / f"symbolic_tra_{expected_code}.svcf",
        svtype="TRA",
        alt="<TRA>",
        chr2=chr2,
        end=end,
    )

    assert expected_code in _codes(_validate(path))


def test_v11_symbolic_tra_requires_dot_svlen(tmp_path):
    path = _write_v11_event(
        tmp_path / "symbolic_tra_bad_svlen.svcf",
        svtype="TRA",
        alt="<TRA>",
        svlen="4900",
    )

    assert "E_TRA_006" in _codes(_validate(path))


def test_v11_symbolic_bnd_remains_invalid(tmp_path):
    path = _write_v11_event(
        tmp_path / "symbolic_bnd.svcf",
        svtype="BND",
        alt="<BND>",
    )

    assert "E_TRA_003" in _codes(_validate(path))


def test_v11_bracket_bnd_with_matching_chr2_end_is_valid(tmp_path):
    path = _write_v11_event(
        tmp_path / "bracket_bnd.svcf",
        svtype="BND",
        alt="N]chr5:5000]",
    )

    assert _validate(path).status() == "PASS"
