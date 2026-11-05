"""Memory-bounded sample-mode merge (issue #194).

The in-memory merge parses every record of every input before merging, so a
large cohort (hundreds of whole-genome samples) can exhaust memory.  Records
are only ever compared inside one (SV type, chromosome) group, and inside
one global TRA group and one global BND group, so the merge can run one
partition at a time:

1. Spool: read each input once, in input order, and append every data line
   (with its input index and line number) to the disk bucket of its
   partition: (SVTYPE, CHROM) for ordinary SVs, one bucket for all TRA and
   one for all BND.  No event objects are kept.
2. Merge: for each bucket, build the events with the same line parser as the
   in-memory path, then run the unchanged classifiers, SVMerger (including
   the #192 split), selection and sample-mode preparation.  The bucket's
   records are written in output order to a run file, each with its output
   sort key, and the bucket's objects are released.
3. Assemble: write the header (with the MERGE_SPLIT definition only if some
   record carries it) and merge the sorted runs by the same output sort key.

Grouping, representative selection, sample synthesis and the SVCF 1.1 record
text are the in-memory path's own code; only loading and output assembly
differ.  The output is byte-identical to the in-memory merge.
"""

from __future__ import annotations

import gc
import heapq
import io
import json
from collections import OrderedDict
from contextlib import ExitStack
from pathlib import Path
import tempfile

from octopusv.merger.multi_sample_writer import MultiSampleWriter
from octopusv.merger.member_breakpoint_range import needs_member_range_headers
from octopusv.merger.sv_merge_writer import warn_legacy_sample_inference
from octopusv.merger.sv_merger import SVMerger
from octopusv.utils.SV_classifier_by_chromosome import SVClassifiedByChromosome
from octopusv.utils.SV_classifier_by_type import SVClassifierByType
from octopusv.utils.svcf_parser import (
    event_from_parts,
    parse_info_field,
    sample_names_from_header,
    split_record_line,
)
from octopusv.utils.svcf_sort import contig_rank, output_sort_key, sort_events_for_output
from octopusv.utils.text_io import open_text_auto

# Sorted runs merged at once; more runs are merged in rounds.
MERGE_FAN_IN = 32
# Bucket files kept open at once while spooling.
MAX_OPEN_BUCKETS = 32
WHOLE_GENOME_TYPES = ("TRA", "BND")


def partition_key(parts):
    """Partition of one split data line.

    Uses the same INFO parser as SVCFEvent, so SVTYPE is read exactly as the
    classifiers read it.  Ordinary SVs are grouped by SVTYPE and CHROM (their
    start chromosome); TRA and BND are matched across chromosome pairs and
    keep one bucket each.  Unknown types get their own bucket, where the
    classifier rejects them as before.
    """
    svtype = str(parse_info_field(parts[7]).get("SVTYPE", ""))
    if svtype in WHOLE_GENOME_TYPES:
        return (svtype, "")
    return (svtype, parts[0])


class _Spool:
    """Append lines to per-partition files with a bounded number of open handles."""

    def __init__(self, directory: Path):
        self.directory = directory
        self.paths: dict[tuple[str, str], Path] = {}
        self._handles: OrderedDict = OrderedDict()

    def append(self, key, text):
        path = self.paths.get(key)
        if path is None:
            path = self.directory / f"bucket_{len(self.paths):05d}.txt"
            self.paths[key] = path
        handle = self._handles.get(key)
        if handle is None:
            if len(self._handles) >= MAX_OPEN_BUCKETS:
                _key, oldest = self._handles.popitem(last=False)
                oldest.close()
            handle = open(path, "a", encoding="utf-8")
            self._handles[key] = handle
        else:
            self._handles.move_to_end(key)
        handle.write(text)

    def close(self):
        while self._handles:
            _key, handle = self._handles.popitem(last=False)
            handle.close()


def _spool_inputs(input_files, directory):
    """One pass over the inputs; returns bucket paths and per-input sample names."""
    spool = _Spool(directory)
    sample_names = []
    try:
        for index, filename in enumerate(input_files):
            names = None
            with open_text_auto(filename) as handle:
                for line_number, line in enumerate(handle, 1):
                    if line.startswith("#CHROM"):
                        names = sample_names_from_header(line, filename)
                        continue
                    if line.startswith("#") or not line.strip():
                        continue
                    parts = split_record_line(line, filename, line_number)
                    if names is None:
                        raise ValueError(
                            f"Malformed SVCF {filename!r}: data record on line "
                            f"{line_number} appears before a #CHROM header."
                        )
                    # The line is stored as read; it is split again with the
                    # same parser when its bucket is merged.
                    spool.append(partition_key(parts), f"{index}\t{line_number}\t{line.rstrip(chr(10))}\n")
            sample_names.append(names)
    finally:
        spool.close()
    return list(spool.paths.values()), sample_names


def _read_bucket(path, input_files, sample_names):
    events = []
    with open(path, encoding="utf-8") as handle:
        for stored in handle:
            index, line_number, line = stored.rstrip("\n").split("\t", 2)
            filename = input_files[int(index)]
            parts = split_record_line(line, filename, int(line_number))
            events.append(event_from_parts(parts, filename, sample_names[int(index)]))
    return events


def _write_run(path, events, writer, rank, contig_order):
    """Records in output order, each line prefixed by its JSON sort key."""
    has_split = False
    range_headers = set()
    with open(path, "w", encoding="utf-8") as handle:
        for event in sort_events_for_output(events, contig_order):
            has_split = has_split or "MERGE_SPLIT" in getattr(event, "info", {})
            range_headers.update(needs_member_range_headers((event,)))
            buffer = io.StringIO()
            writer.write_record(buffer, event)
            key = json.dumps(output_sort_key(event, rank), separators=(",", ":"))
            handle.write(key + "\t" + buffer.getvalue())
    return has_split, range_headers


def _run_lines(handle):
    for stored in handle:
        key, line = stored.split("\t", 1)
        yield json.loads(key), line


def _merge_runs(paths, target, keep_keys):
    with ExitStack() as stack:
        streams = [_run_lines(stack.enter_context(open(path, encoding="utf-8"))) for path in paths]
        for key, line in heapq.merge(*streams, key=lambda item: item[0]):
            if keep_keys:
                target.write(json.dumps(key, separators=(",", ":")) + "\t" + line)
            else:
                target.write(line)


def _merge_runs_bounded(paths, target, directory):
    """Merge sorted runs into target, at most MERGE_FAN_IN files at a time."""
    current = list(paths)
    level = 0
    while len(current) > MERGE_FAN_IN:
        merged = []
        for start in range(0, len(current), MERGE_FAN_IN):
            path = directory / f"merged_{level}_{start:05d}.txt"
            with open(path, "w", encoding="utf-8") as handle:
                _merge_runs(current[start:start + MERGE_FAN_IN], handle, keep_keys=True)
            merged.append(path)
        for old in current:
            old.unlink()
        current = merged
        level += 1
    _merge_runs(current, target, keep_keys=False)


def merge_samples_partitioned(*, input_files, output_file, contigs, name_mapper,
                              prepare_inputs, merge_options, select):
    """Sample-mode merge of input_files into output_file, one partition at a time.

    prepare_inputs(events) prepares freshly parsed sample-mode inputs (the
    CLI's evidence preservation); select(merger) applies the merge strategy.
    Temporary files live next to output_file and are always removed.
    Returns the number of records written.
    """
    filenames = [str(path) for path in input_files]
    contig_order = list(contigs.keys())
    rank = contig_rank(contig_order)
    writer = MultiSampleWriter(name_mapper, input_files=filenames)
    summary = {"legacy_sample_events": 0, "legacy_sample_unresolved_evidence_blocks": 0}
    record_count = 0
    has_split = False
    member_headers = set()

    parent = Path(output_file).parent
    with tempfile.TemporaryDirectory(prefix=".octopusv_merge_", dir=parent) as temp:
        temp = Path(temp)
        buckets, sample_names = _spool_inputs(filenames, temp)
        runs = []
        for number, bucket in enumerate(buckets):
            events = _read_bucket(bucket, filenames, sample_names)
            bucket.unlink()
            prepare_inputs(events)
            by_type = SVClassifierByType(events)
            by_type.classify()
            by_chrom = SVClassifiedByChromosome(by_type.get_classified_events())
            by_chrom.classify()
            merger = SVMerger(by_chrom.get_classified_events(), all_input_files=filenames,
                              mode="sample", **merge_options)
            merger.merge()
            prepared, bucket_summary = merger.prepare_sample_mode_events(select(merger), name_mapper)
            for key in summary:
                summary[key] += bucket_summary.get(key, 0)
            record_count += len(prepared)
            run = temp / f"run_{number:05d}.txt"
            bucket_has_split, bucket_headers = _write_run(run, prepared, writer, rank, contig_order)
            has_split = bucket_has_split or has_split
            member_headers.update(bucket_headers)
            runs.append(run)
            del events, by_type, by_chrom, merger, prepared
            gc.collect()

        with open(output_file, "w") as target:
            writer.write_header(target, contigs, merge_split=has_split,
                                member_ranges=tuple(key for key in ("MERGE_POS_RANGE", "MERGE_END_RANGE")
                                                    if key in member_headers))
            _merge_runs_bounded(runs, target, temp)

    warn_legacy_sample_inference(summary)
    return record_count
