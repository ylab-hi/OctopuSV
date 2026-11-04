"""Split merged candidate groups that hold different structures from one source.

Issue #192.  Core grouping compares every new event with a group's seed only
and does not look at where events come from.  Two records from the same source
can therefore enter one group even though they describe different structures,
and one of them then disappears behind a single representative.

This module runs after core grouping and before representative selection.  It
leaves every group without such a conflict exactly as it is.  A conflicting
group is partitioned so that no part holds two non-equivalent records from one
source, and every other record joins only a part it is equivalent to (or
starts its own); all parts carry the same INFO/MERGE_SPLIT value so the link
is kept.

Equivalence of two records from the same source:

* caller mode (one input file per caller): DEL, DUP and INS are equivalent
  when their lengths differ by at most max(10 bp, 10%).  INV, TRA and BND are
  not changed here: a caller writing one inversion as two halves, or one
  junction as two ends, is the expected case and stays merged.
* sample mode (one input file per sample): the two records must also merge
  directly with each other under the normal matching rule of their type;
  DEL, DUP and INS must in addition pass the length rule.

A length difference above the limit is a reason to keep records apart, not a
proof that they are distinct alleles.
"""

from __future__ import annotations

from collections import defaultdict

from octopusv.utils.source_path import normalize_source_path

LENGTH_RULE_TYPES = {"DEL", "DUP", "INS"}
LENGTH_ABS_BP = 10
LENGTH_REL = 0.10
MODE_CALLER = "caller"
MODE_SAMPLE = "sample"


def event_length(event):
    """Linear length used by the length rule; mirrors should_merge()."""
    svlen = event.info.get("SVLEN", ".")
    try:
        if svlen == "." or not svlen:
            return event.end_pos - event.start_pos + 1
        return abs(int(svlen))
    except (ValueError, TypeError):
        return event.end_pos - event.start_pos + 1


def lengths_equivalent(length_a, length_b):
    longer = max(length_a, length_b)
    return abs(length_a - length_b) <= max(LENGTH_ABS_BP, LENGTH_REL * longer)


def event_svtype(event):
    svtype = getattr(event, "sv_type", None)
    if svtype:
        return svtype
    return str(getattr(event, "info", {}).get("SVTYPE", ""))


def source_key(event):
    original = getattr(event, "_original_source_file", None)
    return normalize_source_path(original if original is not None else event.source_file)


MERGE_SPLIT_HEADER = (
    '##INFO=<ID=MERGE_SPLIT,Number=1,Type=String,Description="Records sharing this '
    'value were one candidate merge group, kept apart because one source reported '
    'records that are not the same structure (DEL/DUP/INS length differs by more '
    'than max(10 bp, 10%); in sample mode also records that do not merge directly). '
    'The value is CHROM:POS:SVTYPE:ID of the group seed">'
)


def needs_split_header(events):
    return any("MERGE_SPLIT" in getattr(event, "info", {}) for event in events)


def split_key(seed):
    """INFO/MERGE_SPLIT value shared by all parts of one split group."""
    return f"{seed.chrom}:{seed.pos}:{event_svtype(seed)}:{seed.sv_id}"


class GroupSplitter:
    """Partition one candidate group; see the module docstring for the rules."""

    def __init__(self, mode, direct_match=None):
        if mode not in (MODE_CALLER, MODE_SAMPLE):
            raise ValueError(f"Unknown merge mode {mode!r}")
        self.mode = mode
        self.direct_match = direct_match

    def equivalent(self, a, b):
        svtype = event_svtype(a)
        if self.mode == MODE_CALLER:
            if svtype not in LENGTH_RULE_TYPES:
                return True
            return lengths_equivalent(event_length(a), event_length(b))
        if self.direct_match is not None:
            try:
                if not (self.direct_match(a, b) and self.direct_match(b, a)):
                    return False
            except (TypeError, ValueError):
                return False
        if svtype in LENGTH_RULE_TYPES:
            return lengths_equivalent(event_length(a), event_length(b))
        return True

    def _distance(self, a, b):
        length_gap = (
            abs(event_length(a) - event_length(b))
            if event_svtype(a) in LENGTH_RULE_TYPES else 0
        )
        position_gap = abs(a.start_pos - b.start_pos) + abs(a.end_pos - b.end_pos)
        return (0 if self.equivalent(a, b) else 1, length_gap, position_gap)

    def _equivalent_parts(self, group, parts, index):
        return [
            number for number, part in enumerate(parts)
            if self.equivalent(group[part[0]], group[index])
        ]

    def split(self, group):
        """Return [group] when there is no conflict, else the list of parts."""
        if len(group) < 2:
            return [group]
        if self.mode == MODE_CALLER and event_svtype(group[0]) not in LENGTH_RULE_TYPES:
            return [group]

        sources = [source_key(event) for event in group]
        by_source = defaultdict(list)
        for index, source in enumerate(sources):
            by_source[source].append(index)

        anchors = set()
        for indices in by_source.values():
            for i_pos, i in enumerate(indices):
                for j in indices[i_pos + 1:]:
                    if not self.equivalent(group[i], group[j]):
                        anchors.update((i, j))
        if not anchors:
            return [group]

        # Anchors are the records involved in a same-source conflict.  Each
        # joins the first part whose first anchor it is equivalent to and that
        # holds no non-equivalent record of its own source.
        parts = []
        for i in sorted(anchors):
            for part in parts:
                if not self.equivalent(group[part[0]], group[i]):
                    continue
                if any(
                    sources[k] == sources[i] and not self.equivalent(group[k], group[i])
                    for k in part
                ):
                    continue
                part.append(i)
                break
            else:
                parts.append([i])

        # Every other record has no same-source conflict anywhere in the group.
        # It joins the nearest part it is equivalent to.  If it is equivalent
        # to none (for example 150 bp between parts of 100 and 200 bp), it is
        # not forced into one: it starts its own part, which later records
        # equivalent to it may join, and keeps the shared MERGE_SPLIT link.
        for i in range(len(group)):
            if i in anchors:
                continue
            candidates = self._equivalent_parts(group, parts, i)
            if candidates:
                best = min(
                    candidates,
                    key=lambda p: (self._distance(group[i], group[parts[p][0]]), p),
                )
                parts[best].append(i)
            else:
                parts.append([i])

        ordered = sorted((sorted(part) for part in parts), key=lambda part: part[0])
        return [[group[i] for i in part] for part in ordered]


def select_with_split(groups, splitter, select):
    """Representatives for groups, splitting conflicting groups first."""
    representatives = []
    for group in groups:
        parts = splitter.split(group) if splitter is not None else [group]
        if len(parts) == 1:
            representatives.append(select(parts[0]))
            continue
        key = split_key(group[0])
        for part in parts:
            representative = select(part)
            representative.info["MERGE_SPLIT"] = key
            representatives.append(representative)
    return representatives
