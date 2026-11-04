#!/usr/bin/env python3
"""Mutation check for merge guard tests.

Each mutation injects one deliberate bug into OctopuSV in a fresh process and
runs the guard tests that are supposed to catch it.  A mutation that leaves
its tests passing means those tests are not real protection.

Run from the repository root:
    python tests/run_guard_mutations.py
"""

from __future__ import annotations

import subprocess
import sys
import textwrap

# name -> (patch code executed before pytest, [test node ids that must fail])
MUTATIONS = {
    "tra_no_reciprocal_frame": (
        """
        import octopusv.merger.TRA_merge_logic as m
        def same_frame_only(sig1, sig2, tol):
            a_local, a_remote, a_rc = sig1
            b_local, b_remote, b_rc = sig2
            return (a_rc == b_rc and m._same_endpoint(a_local, b_local, tol)
                    and m._same_endpoint(a_remote, b_remote, tol))
        m._same_breakend_adjacency = same_frame_only
        """,
        ["tests/test_merge_guard_sample_mode.py::test_same_sample_two_ends_of_one_translocation_stay_one_record"],
    ),
    "tra_orientation_ignored": (
        """
        import octopusv.merger.TRA_merge_logic as m
        m._breakend_orientation_signature = lambda event: None
        import octopusv.merger.tra_merger as t
        t._breakend_orientation_signature = lambda event: None
        """,
        ["tests/test_issue_191_regression.py::test_default_keeps_orientation_consistency"],
    ),
    "del_length_ratio_tightened": (
        """
        import octopusv.merger.sv_merge_logic as m
        m._get_length_ratio_threshold = lambda sv_type: 1.5
        """,
        ["tests/test_merge_guard_sample_mode.py::test_same_sample_two_different_deletions_status_quo_before_192"],
    ),
}

RUNNER = """
import sys
{patch}
import pytest
sys.exit(pytest.main(["-q", "--no-header", "-p", "no:cacheprovider", *sys.argv[1:]]))
"""


def main() -> int:
    bad = 0
    # Control: every guard test must pass without a mutation, otherwise a
    # "caught" result below would mean nothing.
    all_tests = sorted({test for _, tests in MUTATIONS.values() for test in tests})
    control = subprocess.run(
        [sys.executable, "-c", RUNNER.format(patch=""), *all_tests],
        capture_output=True, text=True,
    )
    if control.returncode != 0:
        print(control.stdout[-2000:])
        print("CONTROL FAILED: guard tests do not pass on unmutated code")
        return 1
    print(f"control: {len(all_tests)} guard test(s) pass on unmutated code")
    for name, (patch, tests) in MUTATIONS.items():
        code = RUNNER.format(patch=textwrap.dedent(patch))
        for test in tests:
            proc = subprocess.run(
                [sys.executable, "-c", code, test],
                capture_output=True, text=True,
            )
            # Only a failing assertion counts; a crash in the patch does not.
            caught = proc.returncode == 1 and " failed" in proc.stdout
            print(f"{'CAUGHT ' if caught else 'MISSED '} {name:32} {test}")
            if not caught:
                bad += 1
    print("ALL MUTATIONS CAUGHT" if bad == 0 else f"{bad} MUTATION(S) MISSED")
    return 0 if bad == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
