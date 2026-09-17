from __future__ import annotations

from pathlib import Path

from typer.testing import CliRunner

from octopusv.cli.cli import app
from octopusv.utils.svcf_sample_parser import parse_svcf_sample_block


CALLER_FORMAT = "GT:AD:LN:ST:QV:TY:ID:SC:REF:ALT:CO"
SAMPLE_FORMAT = "GT:AD:UC:UV:LN:ST:QV:TY:ID:SC:REF:ALT:CO"
VCF_SAMPLE_FORMAT = "GT:AD:DP:UC:UV:LN"


def _write_caller_svcf(path: Path, caller: str, records: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)

    lines = [
        "##fileformat=VCFv4.2",
        "##source=OctopuSV",
        "##contig=<ID=chr1,length=1000000>",
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE",
    ]

    for record in records:
        pos = int(record["pos"])
        end = int(record["end"])
        length = end - pos
        record_id = str(record["id"])
        gt = str(record["gt"])
        ad = str(record.get("ad", "5,5"))

        info = ";".join(
            [
                "SVTYPE=DEL",
                f"END={end}",
                f"SVLEN=-{length}",
                "CHR2=chr1",
                "SUPPORT=10",
                "SVMETHOD=caller",
                "RTID=.",
                "AF=.",
                "STRAND=.",
                "RNAMES=.",
                "PRECISE",
            ]
        )
        sample = (
            f"{gt}:{ad}:{length}:.:60:DEL:{record_id}:{caller}:"
            f"N:<DEL>:chr1_{pos}-chr1_{end}"
        )
        lines.append(
            f"chr1\t{pos}\t{record_id}\tN\t<DEL>\t60\tPASS\t"
            f"{info}\t{CALLER_FORMAT}\t{sample}"
        )

    path.write_text("\n".join(lines) + "\n")


def _invoke(runner: CliRunner, args: list[str]):
    result = runner.invoke(app, args)
    assert result.exit_code == 0, (
        f"Command failed: octopusv {' '.join(args)}\n"
        f"stdout/stderr:\n{result.stdout}\n"
        f"exception: {result.exception!r}"
    )
    return result


def _records(path: Path) -> list[list[str]]:
    return [
        line.split("\t")
        for line in path.read_text().splitlines()
        if line and not line.startswith("#")
    ]


def _record_at(path: Path, pos: int) -> list[str]:
    matches = [record for record in _records(path) if int(record[1]) == pos]
    assert len(matches) == 1
    return matches[0]


def _parse_info(info_text: str) -> dict[str, object]:
    parsed: dict[str, object] = {}
    for item in info_text.split(";"):
        if "=" in item:
            key, value = item.split("=", 1)
            parsed[key] = value
        else:
            parsed[item] = True
    return parsed


def test_dirty_checkpoint_caller_to_sample_to_validate_to_vcf(tmp_path):
    """Lock the A/B/C checkpoint across the real producer/consumer chain."""
    runner = CliRunner()

    # Sample A deliberately uses two different physical inputs with the same
    # basename. One caller contributes two nearby records to the common event,
    # including colon-containing Manta IDs and HET/HOM disagreement.
    sample_a_manta = tmp_path / "sampleA" / "manta" / "shared.svcf"
    sample_a_cute = tmp_path / "sampleA" / "cute" / "shared.svcf"
    _write_caller_svcf(
        sample_a_manta,
        "manta",
        [
            {
                "pos": 100,
                "end": 600,
                "id": "MantaDEL:1:0:1:0:0:0",
                "gt": "0/1",
            },
            {
                "pos": 105,
                "end": 605,
                "id": "MantaDEL:2:0:1:0:0:0",
                "gt": "1/1",
            },
            {
                "pos": 5000,
                "end": 5500,
                "id": "MantaDEL:onlyA:0:1:0:0:0",
                "gt": "0/1",
                "ad": "4,6",
            },
        ],
    )
    _write_caller_svcf(
        sample_a_cute,
        "cuteSV",
        [
            {
                "pos": 102,
                "end": 602,
                "id": "cute.1",
                "gt": "0/1",
                "ad": "8,2",
            },
            {
                "pos": 9000,
                "end": 9500,
                "id": "cute.onlyB",
                "gt": "0/1",
            },
        ],
    )

    sample_a_caller = tmp_path / "sampleA.caller.svcf"
    _invoke(
        runner,
        [
            "merge",
            "-i",
            str(sample_a_manta),
            "-i",
            str(sample_a_cute),
            "-o",
            str(sample_a_caller),
            "--mode",
            "caller",
            "--caller-names",
            "manta,cuteSV",
            "--union",
        ],
    )

    common_a = _record_at(sample_a_caller, 102)
    common_a_info = _parse_info(common_a[7])
    assert common_a_info["PRECISE"] is True
    assert "PRECISE=True" not in common_a[7]
    assert common_a_info["SOURCES"] == "manta,manta,cuteSV"
    assert common_a_info["SOURCE_IDS"] == (
        "MantaDEL:1:0:1:0:0:0,"
        "MantaDEL:2:0:1:0:0:0,"
        "cute.1"
    )

    caller_blocks = [
        parse_svcf_sample_block(common_a[8], block)
        for block in common_a[9:]
    ]
    assert [block["SC"] for block in caller_blocks] == [
        "manta",
        "manta",
        "cuteSV",
    ]
    assert [block["ID"] for block in caller_blocks] == [
        "MantaDEL:1:0:1:0:0:0",
        "MantaDEL:2:0:1:0:0:0",
        "cute.1",
    ]

    # Same basename, different directories: --specific must distinguish the
    # physical inputs rather than collapsing them by basename.
    specific_output = tmp_path / "sampleA.manta_specific.svcf"
    _invoke(
        runner,
        [
            "merge",
            "-i",
            str(sample_a_manta),
            "-i",
            str(sample_a_cute),
            "-o",
            str(specific_output),
            "--mode",
            "caller",
            "--caller-names",
            "manta,cuteSV",
            "--specific",
            str(sample_a_manta),
        ],
    )
    specific_records = _records(specific_output)
    assert [int(record[1]) for record in specific_records] == [5000]
    specific_info = _parse_info(specific_records[0][7])
    assert specific_info["SOURCES"] == "manta"
    assert specific_info["SOURCE_IDS"] == "MantaDEL:onlyA:0:1:0:0:0"

    # Sample B has the same common biological event. Its first caller uses a
    # partially missing AD to keep that dirty input shape in the checkpoint.
    sample_b_manta = tmp_path / "sampleB" / "manta.svcf"
    sample_b_cute = tmp_path / "sampleB" / "cute.svcf"
    _write_caller_svcf(
        sample_b_manta,
        "manta",
        [
            {
                "pos": 101,
                "end": 601,
                "id": "MantaDEL:B1:0:1:0:0:0",
                "gt": "0/1",
                "ad": "5,.",
            },
            {
                "pos": 15000,
                "end": 15500,
                "id": "MantaDEL:Bonly:0:1:0:0:0",
                "gt": "0/1",
                "ad": "3,7",
            },
        ],
    )
    _write_caller_svcf(
        sample_b_cute,
        "cuteSV",
        [
            {
                "pos": 103,
                "end": 603,
                "id": "cute.B1",
                "gt": "0/1",
                "ad": "7,3",
            }
        ],
    )

    sample_b_caller = tmp_path / "sampleB.caller.svcf"
    _invoke(
        runner,
        [
            "merge",
            "-i",
            str(sample_b_manta),
            "-i",
            str(sample_b_cute),
            "-o",
            str(sample_b_caller),
            "--mode",
            "caller",
            "--caller-names",
            "manta,cuteSV",
            "--union",
        ],
    )

    population_svcf = tmp_path / "population.svcf"
    _invoke(
        runner,
        [
            "merge",
            "-i",
            str(sample_a_caller),
            "-i",
            str(sample_b_caller),
            "-o",
            str(population_svcf),
            "--mode",
            "sample",
            "--sample-names",
            "sampleA,sampleB",
            "--union",
        ],
    )

    validation = _invoke(runner, ["validate-svcf", str(population_svcf)])
    assert "Mode: sample_multi" in validation.stdout
    assert "Errors: 0" in validation.stdout
    assert "Status: PASS" in validation.stdout

    population_text = population_svcf.read_text()
    assert "##OctopuSV_mode=multi" in population_text
    assert "PRECISE=True" not in population_text

    common_population = _record_at(population_svcf, 102)
    common_population_info = _parse_info(common_population[7])
    assert common_population_info["PRECISE"] is True
    assert common_population_info["SOURCES"] == "sampleA,sampleB"
    assert common_population_info["SOURCE_IDS"] == "cute.1,cute.B1"
    assert common_population[8] == SAMPLE_FORMAT

    sample_a = parse_svcf_sample_block(common_population[8], common_population[9])
    sample_b = parse_svcf_sample_block(common_population[8], common_population[10])

    assert sample_a["GT"] == "1/."
    assert sample_a["AD"] == ".,."
    assert sample_a["UC"] == "2"
    assert sample_a["UV"] == "2"
    assert sample_a["LN"] == "500"

    assert sample_b["GT"] == "0/1"
    assert sample_b["AD"] == ".,."
    assert sample_b["UC"] == "2"
    assert sample_b["UV"] == "2"
    assert sample_b["LN"] == "500"

    population_vcf = tmp_path / "population.vcf"
    _invoke(
        runner,
        [
            "svcf2vcf",
            "-i",
            str(population_svcf),
            "-o",
            str(population_vcf),
        ],
    )

    common_vcf = _record_at(population_vcf, 102)
    assert common_vcf[8] == VCF_SAMPLE_FORMAT
    assert common_vcf[9] == "1/.:.,.:.:2:2:500"
    assert common_vcf[10] == "0/1:.,.:.:2:2:500"

    # A single-sample sample-mode file still carries the multi marker and must
    # remain on the sample conversion path.
    single_sample_svcf = tmp_path / "sampleA.only.svcf"
    _invoke(
        runner,
        [
            "merge",
            "-i",
            str(sample_a_caller),
            "-o",
            str(single_sample_svcf),
            "--mode",
            "sample",
            "--sample-names",
            "sampleA",
            "--union",
        ],
    )
    assert "##OctopuSV_mode=multi" in single_sample_svcf.read_text()

    single_sample_vcf = tmp_path / "sampleA.only.vcf"
    _invoke(
        runner,
        [
            "svcf2vcf",
            "-i",
            str(single_sample_svcf),
            "-o",
            str(single_sample_vcf),
        ],
    )
    single_common = _record_at(single_sample_vcf, 102)
    assert single_common[8] == VCF_SAMPLE_FORMAT
    assert single_common[9] == "1/.:.,.:.:2:2:500"
