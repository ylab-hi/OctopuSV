"""Disk-partitioned caller-mode merge with the unchanged scientific grouping.

Keep only one independently mergeable SV partition in memory at a time. Reuse
sample-mode's input spooling, run sorting, and bounded merge; reuse the exact
caller record serializer and SVMerger. Normal events are keyed by SVTYPE and
start chromosome, while TRA/BND retain their full cross-chromosome partitions.
"""

from __future__ import annotations

import gc
import logging
import tempfile
from pathlib import Path

from octopusv.merger.partitioned_sample_merge import (
    _read_bucket,
    _spool_inputs,
    _write_run,
    _merge_runs_bounded,
)
from octopusv.merger.sv_merge_writer import MergeWriterMixin
from octopusv.merger.sv_merger import SVMerger
from octopusv.utils.SV_classifier_by_chromosome import SVClassifiedByChromosome
from octopusv.utils.SV_classifier_by_type import SVClassifierByType
from octopusv.utils.svcf_sort import contig_rank


class _CallerRecordWriter(MergeWriterMixin):
    """Stream caller-mode records using the normal MergeWriterMixin methods."""

    def __init__(self, filenames, name_mapper):
        self.all_input_files = filenames
        self.name_mapper = name_mapper
        self.source_lookup = self._build_caller_source_lookup()
        self._legacy_caller_event_count = 0
        self._legacy_caller_unresolved_evidence_count = 0

    def write_record(self, handle, event):
        self.write_caller_record(handle, event, self.name_mapper, self.source_lookup)

    def write_header(self, handle, contigs, *, merge_split, member_ranges):
        self._merge_split_header = merge_split
        self._member_range_headers = member_ranges
        self._write_vcf_header(handle, contigs, self.all_input_files)

    def warn_legacy_sources(self):
        if self._legacy_caller_event_count:
            logging.warning(
                "Caller-mode writer used legacy source inference for %d "
                "event(s); %d evidence block(s) could not be mapped and "
                "were omitted.",
                self._legacy_caller_event_count,
                self._legacy_caller_unresolved_evidence_count,
            )


def merge_callers_partitioned(*, input_files, output_file, contigs,
                              name_mapper, merge_options, select):
    """Merge caller-mode SVCFs one independently mergeable partition at a time.

    Temp files live alongside the requested output and are cleaned on success
    or failure. The CLI wraps this in its atomic_output_path as usual.
    """
    filenames = [str(path) for path in input_files]
    contig_order = list(contigs.keys())
    rank = contig_rank(contig_order)
    writer = _CallerRecordWriter(filenames, name_mapper)
    count = 0
    has_split = False
    member_headers = set()

    with tempfile.TemporaryDirectory(prefix=".octopusv_merge_", dir=Path(output_file).parent) as tmp:
        temp = Path(tmp)
        buckets, sample_names = _spool_inputs(filenames, temp)
        runs = []
        for number, bucket in enumerate(buckets):
            events = _read_bucket(bucket, filenames, sample_names)
            bucket.unlink()
            by_type = SVClassifierByType(events)
            by_type.classify()
            by_chrom = SVClassifiedByChromosome(by_type.get_classified_events())
            by_chrom.classify()

            merger = SVMerger(
                by_chrom.get_classified_events(),
                all_input_files=filenames,
                mode="caller",
                **merge_options,
            )
            merger.merge()
            selected = select(merger)
            count += len(selected)

            run = temp / f"run_{number:05d}.txt"
            bucket_split, bucket_ranges = _write_run(run, selected, writer, rank, contig_order)
            has_split = has_split or bucket_split
            member_headers.update(bucket_ranges)
            runs.append(run)
            del events, by_type, by_chrom, merger, selected
            gc.collect()

        with open(output_file, "w", encoding="utf-8") as handle:
            writer.write_header(
                handle, contigs, merge_split=has_split,
                member_ranges=tuple(key for key in ("MERGE_POS_RANGE", "MERGE_END_RANGE")
                                    if key in member_headers),
            )
            _merge_runs_bounded(runs, handle, temp)

    writer.warn_legacy_sources()
    return count
