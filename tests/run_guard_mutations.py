#!/usr/bin/env python3
"""Mutation check for the merge guard tests (run before and during #192).

Each mutation injects one deliberate bug into OctopuSV in a fresh process and
runs all guard tests.  Two requirements:

1. Every mutation must make each test listed for it fail (by assertion, not by
   a crash in the patch).
2. Every guard test must be failed by at least one mutation; a test no
   mutation can break is not real protection.

A control run first checks that all guard tests pass on unmutated code.

Run from the repository root:
    python tests/run_guard_mutations.py
"""

from __future__ import annotations

import subprocess
import sys
import textwrap

GUARD_TESTS = [
    "tests/test_issue_192_merge_split.py",
    "tests/test_merge_guard_vcap_fixtures.py",
    "tests/test_merge_guard_sample_mode.py",
    "tests/test_issue_191_regression.py",
]

F = "tests/test_merge_guard_vcap_fixtures.py::"
S = "tests/test_merge_guard_sample_mode.py::"
T191 = "tests/test_issue_191_regression.py::"

# name -> (patch executed before pytest, tests that must fail)
MUTATIONS = {
    "tra_reciprocal_frame_removed": (
        """
        import octopusv.merger.TRA_merge_logic as m
        def same_frame_only(sig1, sig2, tol):
            a_local, a_remote, a_rc = sig1
            b_local, b_remote, b_rc = sig2
            return (a_rc == b_rc and m._same_endpoint(a_local, b_local, tol)
                    and m._same_endpoint(a_remote, b_remote, tol))
        m._same_breakend_adjacency = same_frame_only
        """,
        [F + "test_both_manta_ends_of_one_translocation_merge",
         F + "test_sv3_translocation_merges_eight_callers_with_same_representative",
         S + "test_same_sample_two_ends_of_one_translocation_stay_one_record"],
    ),
    "tra_orientation_ignored": (
        """
        import octopusv.merger.TRA_merge_logic as m
        import octopusv.merger.tra_merger as t
        m._breakend_orientation_signature = lambda event: None
        t._breakend_orientation_signature = lambda event: None
        """,
        [F + "test_both_manta_ends_of_one_translocation_merge",
         T191 + "test_default_keeps_orientation_consistency"],
    ),
    "tra_strand_default_flipped_off": (
        """
        import inspect
        from octopusv.cli.merge import merge
        inspect.signature(merge).parameters["tra_strand_consistency"].default.default = False
        """,
        [T191 + "test_explicit_flag_is_identical_to_default"],
    ),
    "length_ratio_limit_1_5": (
        """
        import octopusv.merger.sv_merge_logic as m
        m._get_length_ratio_threshold = lambda sv_type: 1.5
        """,
        [F + "test_fixture_output_is_locked[del_chr5_56_112_a]"],
    ),
    "inv_distance_5bp": (
        """
        import octopusv.merger.sv_merge_logic as m
        original = m._get_distance_threshold
        m._get_distance_threshold = lambda t, l: 5 if t == "INV" else original(t, l)
        """,
        [F + "test_two_delly_halves_of_one_inversion_merge",
         F + "test_inv_00000658_status_quo_pending_review"],
    ),
    "dup_length_ratio_exact": (
        """
        import octopusv.merger.sv_merge_logic as m
        original = m._get_length_ratio_threshold
        m._get_length_ratio_threshold = lambda t: 1.00001 if t == "DUP" else original(t)
        """,
        [F + "test_pbsv_split_dup_and_bnd_dup_merge",
         F + "test_large_dup_with_two_records_from_one_caller_merges"],
    ),
    "bnd_delta_widened": (
        """
        from octopusv.merger.bnd_merger import BNDMerger
        original = BNDMerger.__init__
        def wide(self, delta=50, **kwargs):
            original(self, 2000, **kwargs)
        BNDMerger.__init__ = wide
        """,
        [F + "test_manta_single_inverted_junction_stays_two_bnd_records"],
    ),
    "representative_by_caller_rank_only": (
        """
        import functools
        import octopusv.merger.sv_selector as s
        import octopusv.merger.sv_merger as a
        import octopusv.merger.tra_merger as b
        import octopusv.merger.bnd_merger as c
        f = functools.partial(s.select_representative_sv,
                              weights={"support": 0, "quality": 0, "caller": 1})
        a.select_representative_sv = b.select_representative_sv = c.select_representative_sv = f
        """,
        [F + "test_sv3_translocation_merges_eight_callers_with_same_representative",
         F + "test_chr5_56bp_and_112bp_deletions_are_kept_apart"],
    ),
    "evidence_order_reversed": (
        """
        from octopusv.merger.sv_merge_writer import MergeWriterMixin as W
        original = W._prepare_caller_records_exact
        def reversed_records(self, event, name_mapper=None):
            records = original(self, event, name_mapper=name_mapper)
            return None if records is None else list(reversed(records))
        W._prepare_caller_records_exact = reversed_records
        """,
        [F + "test_both_manta_ends_of_one_translocation_merge",
         F + "test_two_delly_halves_of_one_inversion_merge",
         T191 + "test_no_flag_disables_orientation_consistency"],
    ),
    "split_disabled": (
        """
        import octopusv.merger.merge_split as m
        m.GroupSplitter.split = lambda self, group: [group]
        """,
        [F + "test_chr5_56bp_and_112bp_deletions_are_kept_apart",
         F + "test_chr5_second_56_112_site_is_kept_apart",
         S + "test_same_sample_two_different_deletions_are_kept_apart"],
    ),
    "split_length_rule_loosened_to_2_1x": (
        """
        import octopusv.merger.merge_split as m
        m.lengths_equivalent = lambda a, b: max(a, b) <= 2.1 * min(a, b)
        """,
        [F + "test_chr5_56bp_and_112bp_deletions_are_kept_apart"],
    ),
    "split_ignores_direct_match_in_sample_mode": (
        """
        import octopusv.merger.merge_split as m
        original = m.GroupSplitter.__init__
        def init(self, mode, direct_match=None):
            original(self, mode, None)
        m.GroupSplitter.__init__ = init
        m.lengths_equivalent = lambda a, b: True
        """,
        [S + "test_same_sample_two_different_deletions_are_kept_apart"],
    ),
    "split_also_applied_to_tra_in_caller_mode": (
        """
        import octopusv.merger.merge_split as m
        m.LENGTH_RULE_TYPES = {"DEL", "DUP", "INS", "TRA", "INV"}
        m.lengths_equivalent = lambda a, b: a == b
        """,
        [F + "test_both_manta_ends_of_one_translocation_merge"],
    ),
    "split_forces_others_into_nearest_part": (
        """
        import octopusv.merger.merge_split as m
        m.GroupSplitter._equivalent_parts = lambda self, group, parts, i: list(range(len(parts)))
        """,
        ["tests/test_issue_192_merge_split.py::"
         "test_record_between_two_split_structures_is_not_forced_into_either"],
    ),
    "split_never_joins_by_normal_rule": (
        """
        import octopusv.merger.merge_split as m
        m.GroupSplitter._clearly_nearest_part = lambda self, group, parts, i, limit: None
        """,
        ["tests/test_issue_192_merge_split.py::"
         "test_other_caller_joins_the_one_part_it_normally_merges_with"],
    ),
    "split_forces_indistinguishable_record_into_nearest": (
        """
        import octopusv.merger.merge_split as m
        m.GroupSplitter._clearly_nearest_part = lambda self, group, parts, i, limit: min(
            range(limit), key=lambda p: (self._distance(group[i], group[parts[p][0]])[1:], p))
        """,
        ["tests/test_issue_192_merge_split.py::"
         "test_record_between_two_split_structures_is_not_forced_into_either"],
    ),
    "group_compares_with_newest_member": (
        """
        from octopusv.merger.sv_merger import SVMerger
        from octopusv.merger.sv_merge_logic import should_merge
        def add(self, sv_type, chromosome, new_event):
            self._all_merged_cache = None
            events = self.merged_events[sv_type][chromosome]
            groups = self.event_groups[sv_type][chromosome]
            for idx in range(len(events)):
                if should_merge(groups[idx][-1], new_event, self.max_distance,
                                self.max_length_ratio, self.min_jaccard):
                    groups[idx].append(new_event)
                    return
            events.append(new_event)
            groups.append([new_event])
        SVMerger._add_and_merge_sorted_event = add
        """,
        [F + "test_fixture_output_is_locked[del_chr5_56_112_a]"],
    ),
}

RUNNER = """
import sys
{patch}
import pytest
sys.exit(pytest.main(["-q", "--no-header", "-p", "no:cacheprovider", "-rf", *sys.argv[1:]]))
"""


def _run(patch: str):
    proc = subprocess.run(
        [sys.executable, "-c", RUNNER.format(patch=textwrap.dedent(patch)), *GUARD_TESTS],
        capture_output=True, text=True,
    )
    failed = {
        line.split(" ", 2)[1]
        for line in proc.stdout.splitlines()
        if line.startswith("FAILED ")
    }
    crashed = proc.returncode not in (0, 1)
    return proc, failed, crashed


def _collect() -> list[str]:
    proc = subprocess.run(
        [sys.executable, "-m", "pytest", "--collect-only", "-q", "-p", "no:cacheprovider",
         *GUARD_TESTS],
        capture_output=True, text=True,
    )
    return [line for line in proc.stdout.splitlines() if "::" in line]


def main() -> int:
    problems = 0
    proc, failed, crashed = _run("")
    if crashed or failed:
        print(proc.stdout[-3000:])
        print("CONTROL FAILED: guard tests do not pass on unmutated code")
        return 1
    all_tests = _collect()
    print(f"control: {len(all_tests)} guard tests pass on unmutated code\n")

    caught_by = {test: [] for test in all_tests}
    for name, (patch, expected) in MUTATIONS.items():
        proc, failed, crashed = _run(patch)
        if crashed:
            print(f"ERROR   {name}: patch crashed\n{proc.stdout[-1500:]}{proc.stderr[-1500:]}")
            problems += 1
            continue
        missing = [test for test in expected if test not in failed]
        status = "CAUGHT " if not missing else "MISSED "
        print(f"{status} {name:36} {len(failed)} test(s) failed")
        for test in missing:
            print(f"          expected to fail but passed: {test}")
        problems += bool(missing)
        for test in failed:
            caught_by.setdefault(test, []).append(name)

    uncovered = [test for test, names in caught_by.items() if not names]
    print()
    if uncovered:
        print("Guard tests that no mutation breaks:")
        for test in uncovered:
            print(f"  {test}")
        problems += len(uncovered)
    else:
        print(f"Every guard test is broken by at least one mutation ({len(all_tests)} tests).")
    print("ALL MUTATIONS CAUGHT" if problems == 0 else f"{problems} PROBLEM(S)")
    return 0 if problems == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
