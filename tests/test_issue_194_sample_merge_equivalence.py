"""Issue #194: sample-mode merge must give byte-identical output however the
inputs are loaded.

The oracle is the legacy in-memory pipeline, run in-process with exactly the
CLI's defaults, contigs and name mapper: parse every input, classify, merge,
select, write.  `octopusv merge --mode sample` must reproduce its output byte
for byte (only the run-time ##fileDate line is ignored).

The inputs come from real VCaP caller slices (tests/data/merge_guard_vcap):
each caller's records from all fixtures are concatenated, and each "sample"
is a caller-mode merge of a different subset of those callers.  This gives
six samples over ~10 chromosomes with DEL/DUP/INS/INV/TRA/BND, reciprocal TRA
ends, and #192 MERGE_SPLIT sites.  Extra cases cover an
undeclared contig, gzip input, a header-only input, a record with trailing
whitespace, and every selection strategy.
"""

from __future__ import annotations

import gzip
import inspect
from pathlib import Path

import pytest
from typer.testing import CliRunner

from octopusv.cli.cli import app
from octopusv.cli.merge import _preflight_merge_inputs, _preserve_sample_input_evidence, merge
from octopusv.merger.name_mapper import NameMapper
from octopusv.merger.sv_merger import SVMerger
from octopusv.utils.SV_classifier_by_chromosome import SVClassifiedByChromosome
from octopusv.utils.SV_classifier_by_type import SVClassifierByType
from octopusv.utils.svcf_parser import SVCFFileEventCreator


RUNNER = CliRunner()
FIXTURES = Path(__file__).parent / "data" / "merge_guard_vcap"
CALLERS = ["cutesv_supp1", "debreak_supp2", "delly_default", "hitsv_default", "lumpy_default",
           "manta_default", "pbsv_default", "sniffles2_supp2", "svaba_default", "svim_supp1"]
MERGE_OPTIONS = ["tra_delta", "tra_min_overlap_ratio", "tra_strand_consistency",
                 "max_distance", "max_length_ratio", "min_jaccard", "bnd_delta"]


def _cli_default(name):
    default = inspect.signature(merge).parameters[name].default
    return getattr(default, "default", default)


def _caller_files(tmp_path: Path) -> list[Path]:
    """One file per caller: header of its first slice, all distinct records."""
    files = []
    for caller in CALLERS:
        slices = sorted(FIXTURES.glob(f"*/inputs/vcap_{caller}.svcf"))
        assert slices, caller
        header = [line for line in slices[0].read_text().splitlines() if line.startswith("#")]
        records, seen = [], set()
        for path in slices:
            for line in path.read_text().splitlines():
                if line and not line.startswith("#") and line not in seen:
                    seen.add(line)
                    records.append(line)
        out = tmp_path / "callers" / f"{caller}.svcf"
        out.parent.mkdir(exist_ok=True)
        out.write_text("\n".join(header + records) + "\n")
        files.append(out)
    return files


# Each sample is a caller-mode merge of a different subset of real callers,
# as in the real workflow (per-sample caller merge, then sample merge).
SAMPLE_CALLERS = [
    (0, 1, 2, 3, 4, 5, 6, 7, 8, 9),
    (0, 2, 4, 6, 8),
    (1, 3, 5, 7, 9),
    (2, 5, 6, 8),
    (0, 3, 5, 9),
    (6, 7),
]


def _build_samples(tmp_path: Path) -> list[Path]:
    callers = _caller_files(tmp_path)
    samples = []
    for number, subset in enumerate(SAMPLE_CALLERS):
        out = tmp_path / f"sample_{number}.svcf"
        args = ["merge", "--mode", "caller", "--union", "-o", str(out)]
        for index in subset:
            args += ["-i", str(callers[index])]
        result = RUNNER.invoke(app, args)
        assert result.exit_code == 0, result.output
        samples.append(out)
    return samples


def _legacy(inputs: list[Path], selection: dict) -> str:
    """The pre-#194 in-memory sample merge, exactly as the CLI ran it."""
    files = [str(path) for path in inputs]
    mapper = NameMapper(inputs, mode="sample", custom_names=None)
    contigs = _preflight_merge_inputs(
        input_files=inputs, labels=mapper.get_all_display_names(), mode="sample",
        specific=selection.get("specific"), expression=selection.get("expression"),
    )
    creator = SVCFFileEventCreator(files)
    creator.parse()
    _preserve_sample_input_evidence(creator.events)
    by_type = SVClassifierByType(creator.events)
    by_type.classify()
    by_chrom = SVClassifiedByChromosome(by_type.get_classified_events())
    by_chrom.classify()
    merger = SVMerger(by_chrom.get_classified_events(), all_input_files=files, mode="sample",
                      **{name: _cli_default(name) for name in MERGE_OPTIONS})
    merger.merge()
    if "expression" in selection:
        results = merger.get_events_by_expression(selection["expression"])
    elif "intersect" in selection:
        results = merger.get_events_by_source(files, operation="intersection")
    elif "union" in selection:
        results = merger.get_events_by_source(files, operation="union")
    elif "specific" in selection:
        results = merger.get_events_by_source([str(p) for p in selection["specific"]],
                                              operation="specific")
    elif "exact_support" in selection:
        results = merger.get_events_by_exact_support(selection["exact_support"])
    else:
        results = merger.get_events_by_support_range(selection.get("min_support"),
                                                     selection.get("max_support"))
    out = inputs[0].parent / "legacy.svcf"
    merger.write_results(out, results, contigs, "sample", mapper, files)
    return out.read_text()


def _cli(inputs: list[Path], selection: dict) -> str:
    out = inputs[0].parent / "cli.svcf"
    args = ["merge", "--mode", "sample", "-o", str(out)]
    for path in inputs:
        args += ["-i", str(path)]
    for key, value in selection.items():
        flag = "--" + key.replace("_", "-")
        if value is True:
            args.append(flag)
        elif isinstance(value, list):
            for item in value:
                args += [flag, str(item)]
        else:
            args += [flag, str(value)]
    result = RUNNER.invoke(app, args)
    assert result.exit_code == 0, result.output
    assert not list(out.parent.glob(".octopusv_merge_*")), "temporary files left behind"
    validation = RUNNER.invoke(app, ["validate-svcf", "-i", str(out)])
    assert validation.exit_code == 0, validation.output
    return out.read_text()


def _stable(text: str) -> str:
    return "\n".join(line for line in text.splitlines() if not line.startswith("##fileDate="))


def _assert_same(inputs, selection):
    expected = _stable(_legacy(inputs, selection))
    actual = _stable(_cli(inputs, selection))
    assert actual == expected
    return actual


@pytest.mark.parametrize("selection", [
    {"union": True},
    {"intersect": True},
    {"min_support": 2},
    {"max_support": 1},
    {"min_support": 2, "max_support": 4},
    {"exact_support": 3},
    {"expression": "sample_1_svcf AND NOT sample_2_svcf"},
], ids=["union", "intersect", "min2", "max1", "range2_4", "exact3", "expression"])
def test_sample_merge_matches_legacy_for_every_strategy(tmp_path, selection):
    samples = _build_samples(tmp_path)
    text = _assert_same(samples, selection)
    if "union" in selection:
        records = [line for line in text.splitlines() if not line.startswith("#")]
        assert len(records) > 50
        assert any("SVTYPE=TRA" in line for line in records)
        assert any("SVTYPE=BND" in line for line in records)
        assert any("MERGE_SPLIT=" in line for line in records)
        assert len({line.split("\t")[0] for line in records}) >= 8


def test_specific_selection_matches_legacy(tmp_path):
    samples = _build_samples(tmp_path)
    _assert_same(samples, {"specific": samples[5:7]})


def test_undeclared_contig_gzip_empty_input_and_trailing_whitespace(tmp_path):
    samples = _build_samples(tmp_path)

    # Record on a contig no header declares: sorted after declared contigs.
    lines = samples[0].read_text().splitlines()
    donor = next(line for line in lines if not line.startswith("#") and "SVTYPE=DEL" in line)
    fields = donor.split("\t")
    fields[0] = "chrUn_issue194"
    fields[7] = fields[7].replace(f"CHR2={donor.split(chr(9))[0]}", "CHR2=chrUn_issue194")
    lines.append("\t".join(fields))
    # Trailing whitespace after the last sample column.
    trailing = next(i for i, line in enumerate(lines) if not line.startswith("#") and "SVTYPE=INS" in line)
    lines[trailing] = lines[trailing] + "\t"
    samples[0].write_text("\n".join(lines) + "\n")

    # gzip-compressed input.
    gz = samples[1].with_name(samples[1].name + ".gz")
    with gzip.open(gz, "wt") as handle:
        handle.write(samples[1].read_text())
    samples[1].unlink()
    samples[1] = gz

    # Header-only input (a sample without calls).
    header = [line for line in samples[2].read_text().splitlines() if line.startswith("#")]
    empty = tmp_path / "sample_empty.svcf"
    empty.write_text("\n".join(header) + "\n")
    samples.append(empty)

    text = _assert_same(samples, {"union": True})
    assert any(line.startswith("chrUn_issue194\t") for line in text.splitlines())


def test_many_output_runs_use_the_bounded_merge(tmp_path, monkeypatch):
    """Force the multi-level merge of sorted runs (more runs than fan-in)."""
    partitioned = pytest.importorskip("octopusv.merger.partitioned_sample_merge")
    monkeypatch.setattr(partitioned, "MERGE_FAN_IN", 2)
    samples = _build_samples(tmp_path)
    _assert_same(samples, {"union": True})


def _raw_tra_sample(tmp_path: Path, name: str, pos: int, mate: int) -> Path:
    raw = tmp_path / f"{name}.vcf"
    raw.write_text(
        "##fileformat=VCFv4.2\n##source=Sniffles2\n"
        "##contig=<ID=chr1,length=1000000>\n##contig=<ID=chr2,length=1000000>\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n"
        f"chr1\t{pos}\t{name}.tra\tN\tN[chr2:{mate}[\t60\tPASS\tSVTYPE=BND;SUPPORT=8\tGT\t0/1\n"
    )
    svcf = tmp_path / f"{name}.svcf"
    result = RUNNER.invoke(app, ["correct", "-i", str(raw), "-o", str(svcf)])
    assert result.exit_code == 0, result.output
    return svcf


def test_input_order_dependent_tra_grouping_is_kept(tmp_path):
    # A-B and B-C are within the TRA tolerance (2 x 50 bp), A-C is not.  TRA groups compare with the
    # first member in input order, so A,B,C gives {A,B},{C} and C,B,A gives
    # {C,B},{A}.  The partitioned merge must keep the input order.
    a = _raw_tra_sample(tmp_path, "A", 300000, 5000)
    b = _raw_tra_sample(tmp_path, "B", 300060, 5060)
    c = _raw_tra_sample(tmp_path, "C", 300120, 5120)

    forward = _assert_same([a, b, c], {"union": True})
    reverse = _assert_same([c, b, a], {"union": True})

    def groups(text):
        return sorted(
            line.split("\t")[7].split("SOURCE_IDS=")[1].split(";")[0]
            for line in text.splitlines() if not line.startswith("#")
        )

    assert groups(forward) == ["A.tra,B.tra", "C.tra"]
    assert groups(reverse) == ["A.tra", "C.tra,B.tra"]
