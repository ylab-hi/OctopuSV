"""Caller disk-partitioning must equal the historical in-memory merge."""

from __future__ import annotations

import gzip
from pathlib import Path

import pytest

from octopusv.cli.merge import _preflight_merge_inputs
from octopusv.merger.name_mapper import NameMapper
from octopusv.merger.partitioned_caller_merge import merge_callers_partitioned
from octopusv.merger.partitioned_sample_merge import _merge_runs_bounded
from octopusv.merger.sv_merger import SVMerger
from octopusv.utils.SV_classifier_by_chromosome import SVClassifiedByChromosome
from octopusv.utils.SV_classifier_by_type import SVClassifierByType
from octopusv.utils.svcf_parser import SVCFFileEventCreator

FORMAT = "GT:AD:LN:ST:QV:TY:ID:SC:REF:ALT:CO"
OPTIONS = dict(
    tra_delta=50, tra_min_overlap_ratio=0.5, tra_strand_consistency=True,
    max_distance=None, max_length_ratio=None, min_jaccard=0.0, bnd_delta=50,
)


def _record(svtype, chrom, pos, end, record_id, source):
    mate_chrom = "chr2" if svtype in ("TRA", "BND") else chrom
    alt = f"N[{mate_chrom}:{end}[" if svtype in ("TRA", "BND") else f"<{svtype}>"
    svlen = "." if svtype in ("TRA", "BND") else str(end - pos)
    block = (
        f"0/1:5,7:{svlen}:+-:60:{svtype}:{record_id}:{source}:N:{alt}:"
        f"{chrom}_{pos}-{mate_chrom}_{end}"
    )
    return (
        f"{chrom}\t{pos}\t{record_id}\tN\t{alt}\t60\tPASS\t"
        f"SVTYPE={svtype};END={end};SVLEN={svlen};CHR2={mate_chrom};"
        f"SUPPORT=5;SVMETHOD=OctopuSV;RTID=.;AF=.;STRAND=+-;RNAMES=.\t"
        f"{FORMAT}\t{block}\n"
    )


def _fixture(tmp_path):
    files = []
    for index in range(3):
        name = f"caller{index}"
        rows = [
            _record("DEL", "chr1", 100 + index * 3, 160 + index * 3, f"{name}.del", name),
            _record("INS", "chr1", 500 + index * 2, 580 + index * 2, f"{name}.ins", name),
            _record("DUP", "chr3", 800 + index * 3, 980 + index * 3, f"{name}.dup", name),
            _record("INV", "chr4", 1200 + index * 2, 1500 + index * 2, f"{name}.inv", name),
            _record("TRA", "chr1", 3000 + index * 4, 6000 + index * 4, f"{name}.tra", name),
            _record("BND", "chr3", 3200 + index * 3, 6200 + index * 3, f"{name}.bnd", name),
            _record("DEL", "chrUn", 7000 + index * 2, 7100 + index * 2, f"{name}.un", name),
            _record("DEL", "chr2", 20000 + index * 1000, 20090 + index * 1000, f"{name}.far", name),
        ]
        if index == 0:
            # Two distinct records from the same caller in one candidate
            # group exercise #192 and MERGE_SPLIT header management.
            rows.insert(1, _record("DEL", "chr1", 103, 166, f"{name}.del.extra", name))
        header = (
            "##fileformat=VCFv4.2\n##SVCFVersion=1.1\n##OctopuSV_mode=caller\n"
            "##contig=<ID=chr1,length=10000000>\n"
            "##contig=<ID=chr2,length=10000000>\n"
            "##contig=<ID=chr3,length=10000000>\n"
            "##contig=<ID=chr4,length=10000000>\n"
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n"
        )
        f = tmp_path / f"{name}.svcf"
        f.write_text(header + "".join(rows))
        files.append(f)
    return files


def _select(merger, files, strategy):
    if strategy == "union":
        return merger.get_events_by_source(files, operation="union")
    if strategy == "intersect":
        return merger.get_events_by_source(files, operation="intersection")
    if strategy == "specific":
        return merger.get_events_by_source([files[0]], operation="specific")
    if strategy == "min2":
        return merger.get_events_by_support_range(min_support=2)
    if strategy == "range":
        return merger.get_events_by_support_range(min_support=1, max_support=2)
    if strategy == "exact2":
        return merger.get_events_by_exact_support(2)
    if strategy == "expression":
        return merger.get_events_by_expression("caller0_svcf AND NOT caller1_svcf")
    raise AssertionError(strategy)


def _full_merge(files, contigs, mapper, selection, out):
    creator = SVCFFileEventCreator(files)
    creator.parse()
    by_type = SVClassifierByType(creator.events)
    by_type.classify()
    by_chrom = SVClassifiedByChromosome(by_type.get_classified_events())
    by_chrom.classify()
    merger = SVMerger(by_chrom.get_classified_events(), all_input_files=files, mode="caller", **OPTIONS)
    merger.merge()
    result = _select(merger, files, selection)
    merger.write_results(out, result, contigs, "caller", mapper, files)
    return len(result)


def _without_date(path):
    return "".join(line for line in path.read_text().splitlines(keepends=True)
                   if not line.startswith("##fileDate="))


@pytest.mark.parametrize("selection", ["union", "intersect", "specific", "min2", "range", "exact2", "expression"])
@pytest.mark.parametrize("reverse", [False, True])
def test_matches_full_caller_merge_on_all_selection_modes(tmp_path, selection, reverse):
    paths = _fixture(tmp_path)
    if reverse:
        paths.reverse()
    files = [str(p) for p in paths]
    mapper = NameMapper(files, mode="caller")
    contigs = _preflight_merge_inputs(input_files=paths, labels=mapper.get_all_display_names(),
                                       mode="caller", specific=None, expression=None)
    full = tmp_path / "legacy.svcf"
    partitioned = tmp_path / "partitioned.svcf"
    count_full = _full_merge(files, contigs, mapper, selection, full)
    count_partitioned = merge_callers_partitioned(
        input_files=files, output_file=partitioned, contigs=contigs,
        name_mapper=mapper, merge_options=OPTIONS,
        select=lambda merger: _select(merger, files, selection),
    )
    assert count_partitioned == count_full
    assert _without_date(partitioned) == _without_date(full)
    assert not list(tmp_path.glob(".octopusv_merge_*"))


def test_gzip_symlinks_custom_labels_header_only_and_many_run_rounds(tmp_path, monkeypatch):
    from octopusv.merger import partitioned_sample_merge as shared

    files = _fixture(tmp_path)
    gz = tmp_path / "caller1.svcf.gz"
    with gzip.open(gz, "wt") as handle:
        handle.write(files[1].read_text())
    files[1] = gz
    link = tmp_path / "caller2.link.svcf"
    link.symlink_to(files[2])
    files[2] = link
    # A source with no records must still be present in source-name mapping.
    empty = tmp_path / "empty.svcf"
    empty.write_text("\n".join(line for line in files[0].read_text().splitlines() if line.startswith("#")) + "\n")
    files.append(empty)

    # Force multi-round external merge without changing any source metadata.
    monkeypatch.setattr(shared, "MERGE_FAN_IN", 2)
    names = ["A", "B", "C", "EMPTY"]
    paths = [str(p) for p in files]
    mapper = NameMapper(paths, mode="caller", custom_names=names)
    contigs = _preflight_merge_inputs(input_files=files, labels=names, mode="caller",
                                       specific=None, expression=None)
    expected = tmp_path / "legacy.svcf"
    actual = tmp_path / "partitioned.svcf"
    _full_merge(paths, contigs, mapper, "union", expected)
    merge_callers_partitioned(input_files=paths, output_file=actual,
                              contigs=contigs, name_mapper=mapper,
                              merge_options=OPTIONS, select=lambda merger: _select(merger, paths, "union"))
    assert _without_date(actual) == _without_date(expected)
