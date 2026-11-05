"""`correct` must not pass on TRA records whose remote breakpoint is unusable.

Seen in DeBreak output: about half of its TRA records repeat POS as END on
another chromosome (END copied from the local side; the real remote position
is not reported).  Some of these END values are even beyond the length of the
remote chromosome.  Such records cannot describe a translocation (SVCF 1.1
TRA needs two known breakpoints), so `correct` skips them, counts them in the
output header and the summary, and lists every one in a side report.  It does
not guess a coordinate and does not drop anything silently.

Rules:
* any caller: a TRA whose remote position is outside the CHR2 contig length
  -> skipped;
* DeBreak only: symbolic <TRA> on another chromosome with END == POS ->
  skipped (DeBreak's known defect: the local value copied into END);
* everything else (normal DeBreak TRA, bracket-form TRA, same-chromosome
  records, and END == POS from any other caller) is unchanged.
"""

from __future__ import annotations

from pathlib import Path

from typer.testing import CliRunner

from octopusv.cli.cli import app


RUNNER = CliRunner()

CONTIGS = {"chr1": 248956422, "chr8": 145138636, "chr16": 90338345,
           "chr17": 83257441, "chr20": 64444167, "chr21": 46709983}


def _header(source="DeBreak") -> str:
    return (
        "##fileformat=VCFv4.2\n"
        f"##source={source}\n"
        + "".join(f"##contig=<ID={c},length={n}>\n" for c, n in CONTIGS.items())
        + '##ALT=<ID=TRA,Description="Translocation">\n'
        + '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
        + "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tVCaP\n"
    )


def _tra(chrom, pos, record_id, chr2, end, alt="<TRA>", method="DeBreak"):
    return (
        f"{chrom}\t{pos}\t{record_id}\tN\t{alt}\t.\tPASS\t"
        f"END={end};CHR2={chr2};SVLEN=0;SUPPREAD=7;MAPQ=60;SVMETHOD={method};SVTYPE=TRA\tGT\t0/1\n"
    )


def _del(pos, record_id, length=500):
    return (
        f"chr1\t{pos}\t{record_id}\tN\t<DEL>\t.\tPASS\t"
        f"END={pos + length};SVLEN=-{length};SUPPREAD=7;SVMETHOD=DeBreak;SVTYPE=DEL\tGT\t0/1\n"
    )


def _correct(tmp_path, records, source="DeBreak"):
    raw = tmp_path / "debreak.vcf"
    out = tmp_path / "debreak.svcf"
    raw.write_text(_header(source) + "".join(records))
    result = RUNNER.invoke(app, ["correct", "-i", str(raw), "-o", str(out)])
    assert result.exit_code == 0, result.output
    return result, out


def _ids(path: Path) -> set[str]:
    ids = set()
    for line in path.read_text().splitlines():
        if line and not line.startswith("#"):
            ids.add(line.split("\t")[2])
    return ids


def _report(out: Path) -> dict[str, dict[str, str]]:
    report = Path(f"{out}.unresolved_tra.tsv")
    assert report.exists()
    lines = report.read_text().splitlines()
    header = lines[0].split("\t")
    return {row[header.index("ID")]: dict(zip(header, row))
            for row in (line.split("\t") for line in lines[1:])}


RECORDS = [
    # DB7633: END copied from POS and beyond chr17 (83.3 Mb).
    _tra("chr16", 84992413, "DB7633", "chr17", 84992413),
    # DB8415: END copied from POS and beyond chr21 (46.7 Mb).
    _tra("chr17", 75792307, "DB8415", "chr21", 75792307),
    # DB4: END copied from POS, inside chr20, still not a remote coordinate.
    _tra("chr1", 205659, "DB4", "chr20", 205659),
    # DB12 and DB239: normal DeBreak TRA with a real remote position.
    _tra("chr1", 876904, "DB12", "chr20", 29368734),
    _tra("chr1", 8723768, "DB239", "chr8", 144297253),
    # Different END, but beyond chr21: invalid in any form.
    _tra("chr1", 1000000, "BEYOND", "chr21", 50000000),
    # Bracket-form TRA: the remote position is written in ALT, so END == POS
    # here is a coordinate, not a copy.
    _tra("chr1", 3000000, "BRACKET", "chr20", 3000000, alt="N[chr20:3000000["),
    _del(5000000, "DEL1"),
]


def test_unusable_remote_breakpoints_are_skipped_and_reported(tmp_path):
    result, out = _correct(tmp_path, RECORDS)

    assert _ids(out) == {"DB12", "DB239", "BRACKET", "DEL1"}
    validation = RUNNER.invoke(app, ["validate-svcf", "-i", str(out)])
    assert validation.exit_code == 0, validation.output

    report = _report(out)
    assert set(report) == {"DB7633", "DB8415", "DB4", "BEYOND"}
    assert report["DB7633"]["reason"] == "remote_beyond_contig"
    assert report["DB7633"]["CHR2_length"] == str(CONTIGS["chr17"])
    assert report["BEYOND"]["reason"] == "remote_beyond_contig"
    assert report["DB4"]["reason"] == "remote_not_reported"
    assert report["DB4"]["CHR2"] == "chr20" and report["DB4"]["END"] == "205659"

    header = [line for line in out.read_text().splitlines() if line.startswith("##")]
    assert "##OctopuSV_skipped_unresolved_tra=4" in header
    assert "4" in result.output and "unresolved_tra.tsv" in result.output


def test_no_report_when_every_tra_is_usable(tmp_path):
    _result, out = _correct(tmp_path, [_tra("chr1", 876904, "DB12", "chr20", 29368734),
                                       _del(5000000, "DEL1")])

    assert _ids(out) == {"DB12", "DEL1"}
    assert not Path(f"{out}.unresolved_tra.tsv").exists()
    assert "OctopuSV_skipped_unresolved_tra" not in out.read_text()


def test_stale_report_from_an_earlier_run_is_removed(tmp_path):
    _correct(tmp_path, RECORDS)
    _result, out = _correct(tmp_path, [_tra("chr1", 876904, "DB12", "chr20", 29368734)])
    assert not Path(f"{out}.unresolved_tra.tsv").exists()


def test_same_chromosome_record_with_end_equal_pos_is_not_touched(tmp_path):
    # Only END == POS on another chromosome means a copied value.
    _result, out = _correct(tmp_path, [_tra("chr1", 876904, "DB12", "chr20", 29368734),
                                       _tra("chr1", 777777, "SAMECHR", "chr1", 777777)])
    assert "SAMECHR" in _ids(out)
    assert not Path(f"{out}.unresolved_tra.tsv").exists()


def test_without_contig_lengths_only_the_copied_end_rule_applies(tmp_path):
    raw = tmp_path / "debreak.vcf"
    out = tmp_path / "debreak.svcf"
    header = "".join(line + "\n" for line in _header().splitlines()
                     if not line.startswith("##contig"))
    raw.write_text(header + _tra("chr16", 84992413, "DB7633", "chr17", 84992413)
                   + _tra("chr1", 1000000, "BEYOND", "chr21", 50000000)
                   + _tra("chr1", 876904, "DB12", "chr20", 29368734))
    result = RUNNER.invoke(app, ["correct", "-i", str(raw), "-o", str(out)])
    assert result.exit_code == 0, result.output

    # No lengths: BEYOND cannot be judged and is kept; DB7633 is still a copy.
    assert _ids(out) == {"BEYOND", "DB12"}
    assert _report(out)["DB7633"]["reason"] == "remote_not_reported"


def test_other_caller_with_equal_numbers_on_two_chromosomes_is_kept(tmp_path):
    # Another caller reporting chr1:205659 -> chr20:205659 is a legal TRA
    # (two valid breakpoints); only DeBreak's copied END is a known defect.
    # Out-of-range remote positions are still skipped for any caller.
    _result, out = _correct(
        tmp_path,
        [_tra("chr1", 205659, "OTHER_EQUAL", "chr20", 205659, method="OtherCaller"),
         _tra("chr16", 84992413, "OTHER_BEYOND", "chr17", 84992413, method="OtherCaller")],
        source="OtherCaller",
    )

    assert _ids(out) == {"OTHER_EQUAL"}
    assert _report(out)["OTHER_BEYOND"]["reason"] == "remote_beyond_contig"


def test_debreak_is_recognized_by_svmethod_without_source_header(tmp_path):
    _result, out = _correct(tmp_path, [_tra("chr1", 205659, "DB4", "chr20", 205659),
                                       _tra("chr1", 876904, "DB12", "chr20", 29368734)],
                            source="unknown")
    assert _ids(out) == {"DB12"}
    assert _report(out)["DB4"]["reason"] == "remote_not_reported"
