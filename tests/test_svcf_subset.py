from pathlib import Path

from octopusv.subset.svcf_subset import (
    SVCFSubset,
    SubsetConfig,
    split_positional_info_list,
)


FORMAT_FIELD = "GT:AD:LN:ST:QV:TY:ID:SC:REF:ALT:CO"


def _evidence(source_id, caller, pos):
    return (
        f"0/1:5,5:50:.:60:INS:{source_id}:{caller}:N:<INS>:"
        f"chr1_{pos}-chr1_{pos + 50}"
    )


def _write_caller_mode_svcf(
    path: Path,
    *,
    sources: str,
    source_ids: str,
    evidence_blocks: list[str],
):
    path.write_text(
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n"
        "chr1\t100\tmerged.1\tN\t<INS>\t60\tPASS\t"
        "SVTYPE=INS;END=150;SVLEN=50;CHR2=chr1;SUPPORT=10;"
        "SVMETHOD=OctopuSV;RTID=.;AF=.;STRAND=.;RNAMES=.;"
        f"SOURCES={sources};SOURCE_IDS={source_ids}\t"
        f"{FORMAT_FIELD}\t"
        + "\t".join(evidence_blocks)
        + "\n"
    )


def _data_fields(path: Path):
    line = next(
        line
        for line in path.read_text().splitlines()
        if line and not line.startswith("#")
    )
    return line.split("\t")


def test_positional_info_list_preserves_missing_slots():
    assert split_positional_info_list("a,.,c") == ["a", ".", "c"]
    assert split_positional_info_list(".") == ["."]
    assert split_positional_info_list("a,,c") == ["a", ".", "c"]


def test_caller_subset_preserves_other_source_ids_when_one_id_is_missing(tmp_path):
    input_file = tmp_path / "input.svcf"
    output_file = tmp_path / "output.svcf"

    _write_caller_mode_svcf(
        input_file,
        sources="cuteSV,sniffles,svim",
        source_ids="cute.1,.,svim.1",
        evidence_blocks=[
            _evidence("cute.1", "cuteSV", 100),
            _evidence(".", "sniffles", 101),
            _evidence("svim.1", "svim", 102),
        ],
    )

    subset = SVCFSubset(
        SubsetConfig(
            input_file=str(input_file),
            output_file=str(output_file),
            mode="caller",
            selected_callers=["cuteSV", "sniffles", "svim"],
        )
    )
    subset.run()

    fields = _data_fields(output_file)

    assert "SOURCE_IDS=cute.1,.,svim.1" in fields[7]
    assert len(fields[9:]) == 3


def test_caller_statistics_count_unique_sources_not_evidence_blocks(tmp_path):
    input_file = tmp_path / "input.svcf"
    output_file = tmp_path / "output.svcf"

    _write_caller_mode_svcf(
        input_file,
        sources="sniffles,sniffles,pbsv",
        source_ids="sniffles.1,sniffles.2,pbsv.1",
        evidence_blocks=[
            _evidence("sniffles.1", "sniffles", 100),
            _evidence("sniffles.2", "sniffles", 101),
            _evidence("pbsv.1", "pbsv", 102),
        ],
    )

    subset = SVCFSubset(
        SubsetConfig(
            input_file=str(input_file),
            output_file=str(output_file),
            mode="caller",
            selected_callers=["sniffles", "pbsv"],
        )
    )
    subset.run()

    assert subset.available_callers_counter["sniffles"] == 1
    assert subset.available_callers_counter["pbsv"] == 1
    assert subset.retained_callers_counter["sniffles"] == 1
    assert subset.retained_callers_counter["pbsv"] == 1

    fields = _data_fields(output_file)
    assert "SOURCES=sniffles,sniffles,pbsv" in fields[7]
    assert "SOURCE_IDS=sniffles.1,sniffles.2,pbsv.1" in fields[7]
    assert len(fields[9:]) == 3
