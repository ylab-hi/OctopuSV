"""SVCF validator: the OctopuSV SVCF contract checker.

This is NOT a generic VCF lint. It checks whether a file conforms to the
SVCF contract produced by OctopuSV across the three legal modes:

  single-caller   : no ##OctopuSV_mode marker, INFO has no SOURCES,
                    exactly one sample/caller column per row.
  caller-merge    : no ##OctopuSV_mode marker, INFO has SOURCES on every row,
                    number of trailing caller columns == len(SOURCES).
                    Row column count can vary across records by design.
  sample/multi    : has ##OctopuSV_mode=multi, #CHROM lists N sample/caller names,
                    every row has exactly N trailing columns.

Output answers one practical question:
  "Can this SVCF safely enter OctopuSV / downstream analysis / agent workflows,
   and if not, what is broken?"

It emits a human summary by default and structured JSON for --json.
Exit codes:
  0 = pass / pass-with-warnings
  1 = validation failed
  2 = unreadable file
"""

from __future__ import annotations

import json
import os
import re
from dataclasses import dataclass


# OctopuSV's official SVCF FORMAT contract. Order matters.
EXPECTED_FORMAT = "GT:AD:LN:ST:QV:TY:ID:SC:REF:ALT:CO"

# The 9 fixed leading columns of any SVCF/VCF record.
CORE_COLUMNS = [
    "#CHROM",
    "POS",
    "ID",
    "REF",
    "ALT",
    "QUAL",
    "FILTER",
    "INFO",
    "FORMAT",
]

LEGAL_SVTYPES = {"DEL", "DUP", "INV", "INS", "TRA", "BND"}

# INFO keys that OctopuSV writes in current SVCF.
# Value may be ".", but the key should be present.
# SOURCES is mode-dependent and checked separately.
REQUIRED_INFO_KEYS = {
    "SVTYPE",
    "END",
    "SVLEN",
    "CHR2",
    "SUPPORT",
    "SVMETHOD",
    "RTID",
    "AF",
    "STRAND",
    "RNAMES",
}

# Mode marker line written by the SVCF writer for sample/multi mode.
MODE_MULTI_MARKER = "##OctopuSV_mode=multi"

# BND/TRA ALT forms:
#   t[p[
#   t]p]
#   [p[t
#   ]p]t
# where p == chrom:pos.
BND_ALT_RE = re.compile(
    r"([ACGTNacgtn]*)([\[\]])([^:\[\]]+):(\d+)([\[\]])([ACGTNacgtn]*)"
)


@dataclass
class Issue:
    """A single validation finding."""

    level: str
    code: str
    message: str
    record_id: str | None = None
    line_no: int | None = None
    blocking: bool = False

    def to_dict(self) -> dict:
        return {
            "level": self.level,
            "code": self.code,
            "message": self.message,
            "record_id": self.record_id,
            "line_no": self.line_no,
            "blocking_for_downstream": self.blocking,
        }


def _parse_info(info_str: str) -> dict:
    """Parse the INFO column into a dict. Flags without '=' map to True."""
    info = {}
    for item in info_str.split(";"):
        if not item:
            continue
        if "=" in item:
            key, value = item.split("=", 1)
            info[key] = value
        else:
            info[item] = True
    return info


def _has_sources(info: dict) -> bool:
    """Return True if INFO has a non-empty SOURCES value."""
    value = info.get("SOURCES")
    return value not in (None, "", ".", True)


def _parse_sources(info: dict) -> list[str]:
    """Parse SOURCES=caller1,caller2 into a clean caller list."""
    if not _has_sources(info):
        return []
    return [caller for caller in str(info["SOURCES"]).split(",") if caller]


def _is_positive_int(value: object) -> bool:
    return str(value).isdigit() and int(str(value)) > 0


def _is_nonnegative_int(value: object) -> bool:
    return str(value).isdigit() and int(str(value)) >= 0


def _parse_bnd_alt(alt: str) -> tuple[str, int] | None:
    """Parse a BND/TRA ALT bracket into (mate_chrom, mate_pos).

    Returns None if ALT is not a valid full BND bracket form.

    We intentionally parse only the ALT column. This avoids the real SVCF trap
    where FORMAT sample/caller values may contain ':' inside caller IDs.
    """
    match = BND_ALT_RE.fullmatch(alt)
    if not match:
        return None

    _, _, mate_chrom, mate_pos, _, _ = match.groups()
    try:
        return mate_chrom, int(mate_pos)
    except ValueError:
        return None


class SVCFValidator:
    """Validate one SVCF file against the OctopuSV SVCF contract."""

    def __init__(self, path: str, strict_co: bool = False):
        self.path = path
        self.strict_co = strict_co

        self.issues: list[Issue] = []
        self.records_total = 0
        self.mode: str | None = None
        self.declared_samples: list[str] = []
        self.unreadable = False

    # ------------------------------------------------------------------
    # Issue helpers
    # ------------------------------------------------------------------

    def _err(
        self,
        code: str,
        msg: str,
        record_id: str | None = None,
        line_no: int | None = None,
        blocking: bool = False,
    ) -> None:
        self.issues.append(
            Issue(
                level="error",
                code=code,
                message=msg,
                record_id=record_id,
                line_no=line_no,
                blocking=blocking,
            )
        )

    def _warn(
        self,
        code: str,
        msg: str,
        record_id: str | None = None,
        line_no: int | None = None,
    ) -> None:
        self.issues.append(
            Issue(
                level="warning",
                code=code,
                message=msg,
                record_id=record_id,
                line_no=line_no,
                blocking=False,
            )
        )

    # ------------------------------------------------------------------
    # Main entry
    # ------------------------------------------------------------------

    def validate(self) -> None:
        """Run all checks and populate self.issues."""
        if not os.path.exists(self.path) or os.path.getsize(self.path) == 0:
            self.unreadable = True
            self._err("E_FILE_001", f"File not found or empty: {self.path}")
            return

        try:
            with open(self.path, encoding="utf-8-sig") as handle:
                lines = handle.readlines()
        except (OSError, UnicodeDecodeError) as exc:
            self.unreadable = True
            self._err("E_FILE_002", f"Cannot read file: {exc}")
            return

        header_lines = [line.rstrip("\n") for line in lines if line.startswith("##")]
        chrom_line = next(
            (line.rstrip("\n") for line in lines if line.startswith("#CHROM")),
            None,
        )
        data_lines = [
            (i + 1, line.rstrip("\n"))
            for i, line in enumerate(lines)
            if line.strip() and not line.startswith("#")
        ]

        self._check_header(chrom_line)
        if chrom_line is None:
            return

        self._detect_mode(header_lines, chrom_line, data_lines)

        for line_no, line in data_lines:
            self.records_total += 1
            self._check_record(line, line_no)

    # ------------------------------------------------------------------
    # Header / mode
    # ------------------------------------------------------------------

    def _check_header(self, chrom_line: str | None) -> None:
        """Verify that #CHROM exists and has the correct 9 leading columns."""
        if chrom_line is None:
            self._err(
                "E_HDR_001",
                "Missing #CHROM header line.",
                blocking=True,
            )
            return

        columns = chrom_line.split("\t")
        if columns[:9] != CORE_COLUMNS:
            self._err(
                "E_HDR_002",
                f"#CHROM leading columns must be {CORE_COLUMNS}, got {columns[:9]}.",
                blocking=True,
            )

    def _detect_mode(
        self,
        header_lines: list[str],
        chrom_line: str,
        data_lines: list[tuple[int, str]],
    ) -> None:
        """Classify the whole file into one legal SVCF mode.

        Official contract:
          marker present                    -> sample_multi
          no marker + all rows have SOURCES -> caller_merge
          no marker + no row has SOURCES    -> single
          no marker + mixed SOURCES status  -> invalid mixed mode
        """
        columns = chrom_line.split("\t")

        if any(line.strip() == MODE_MULTI_MARKER for line in header_lines):
            self.mode = "sample_multi"
            self.declared_samples = columns[9:]

            if not self.declared_samples:
                self._err(
                    "E_MODE_001",
                    "sample/multi mode but #CHROM declares no sample columns.",
                    blocking=True,
                )
            return

        source_flags = []
        malformed_for_mode = 0

        for _, line in data_lines:
            parts = line.split("\t")
            if len(parts) < 8:
                malformed_for_mode += 1
                continue
            info = _parse_info(parts[7])
            source_flags.append(_has_sources(info))

        if not source_flags:
            self.mode = "single"
            if malformed_for_mode == 0:
                self._warn(
                    "W_MODE_001",
                    "No data records found; treating file as single-caller mode.",
                )
            return

        has_sources_count = sum(source_flags)
        no_sources_count = len(source_flags) - has_sources_count

        if has_sources_count and no_sources_count:
            self.mode = "mixed"
            self._err(
                "E_MODE_002",
                "Invalid mixed SVCF mode: some records have SOURCES and some do not. "
                "A no-marker SVCF must be either all single-caller or all caller-merge.",
                blocking=True,
            )
        elif has_sources_count:
            self.mode = "caller_merge"
        else:
            self.mode = "single"

    # ------------------------------------------------------------------
    # Per-record checks
    # ------------------------------------------------------------------

    def _check_record(self, line: str, line_no: int) -> None:
        parts = line.split("\t")
        if len(parts) < 10:
            self._err(
                "E_REC_001",
                f"Record has fewer than 10 columns ({len(parts)}).",
                line_no=line_no,
                blocking=True,
            )
            return

        chrom, pos, sv_id = parts[0], parts[1], parts[2]
        alt, info_str, fmt = parts[4], parts[7], parts[8]
        sample_cols = parts[9:]
        info = _parse_info(info_str)

        self._check_core_fields(chrom, pos, sv_id, line_no)
        self._check_format(fmt, sv_id, line_no)
        self._check_info_keys(info, sv_id, line_no)
        self._check_support(info, sv_id, line_no)
        self._check_sources_by_mode(info, sv_id, line_no)
        self._check_column_count(info, sample_cols, sv_id, line_no)
        self._check_svtype_and_coords(info, alt, pos, sv_id, line_no)
        self._check_co(sample_cols, sv_id, line_no)

    def _check_core_fields(
        self,
        chrom: str,
        pos: str,
        sv_id: str,
        line_no: int,
    ) -> None:
        if not chrom:
            self._err(
                "E_CHROM_001",
                "CHROM is empty.",
                sv_id,
                line_no,
                blocking=True,
            )

        if not _is_positive_int(pos):
            self._err(
                "E_POS_001",
                f"POS must be a positive integer, got '{pos}'.",
                sv_id,
                line_no,
                blocking=True,
            )

    def _check_format(self, fmt: str, sv_id: str, line_no: int) -> None:
        """FORMAT must equal the official SVCF contract exactly."""
        if fmt != EXPECTED_FORMAT:
            self._err(
                "E_FMT_001",
                f"FORMAT must be '{EXPECTED_FORMAT}', got '{fmt}'.",
                sv_id,
                line_no,
                blocking=True,
            )

    def _check_info_keys(self, info: dict, sv_id: str, line_no: int) -> None:
        """Required INFO keys must be present. Values may be '.'."""
        missing = REQUIRED_INFO_KEYS - set(info.keys())
        if missing:
            blocking = bool({"SVTYPE", "END", "CHR2"} & missing)
            self._err(
                "E_INFO_001",
                f"Missing required INFO field(s): {sorted(missing)}.",
                sv_id,
                line_no,
                blocking=blocking,
            )

    def _check_support(self, info: dict, sv_id: str, line_no: int) -> None:
        """SUPPORT is read support, not caller count.

        OctopuSV SVCF allows SUPPORT='.' when record-level read support is
        unavailable or not applicable.

        Contract:
          - missing SUPPORT key is handled by _check_info_keys;
          - SUPPORT='.' is valid;
          - numeric SUPPORT must be a non-negative integer;
          - any other value is invalid.
        """
        support = info.get("SUPPORT")

        if support is None:
            return

        if support == ".":
            return

        if support is True or support == "" or not _is_nonnegative_int(support):
            self._err(
                "E_SUPPORT_001",
                f"SUPPORT must be '.' or a non-negative integer read-support value, "
                f"got '{support}'.",
                sv_id,
                line_no,
                blocking=False,
            )

    def _check_sources_by_mode(
        self,
        info: dict,
        sv_id: str,
        line_no: int,
    ) -> None:
        """Check SOURCES presence according to official SVCF mode."""
        if self.mode == "single":
            if _has_sources(info):
                self._err(
                    "E_SRC_001",
                    "single-caller mode must not contain SOURCES.",
                    sv_id,
                    line_no,
                    blocking=True,
                )

        elif self.mode in {"caller_merge", "sample_multi"}:
            if not _has_sources(info):
                self._err(
                    "E_SRC_002",
                    f"{self.mode} mode requires non-empty SOURCES.",
                    sv_id,
                    line_no,
                    blocking=True,
                )

    def _check_column_count(
        self,
        info: dict,
        sample_cols: list[str],
        sv_id: str,
        line_no: int,
    ) -> None:
        """Mode-aware sample/caller column-count contract."""
        n_cols = len(sample_cols)

        if self.mode == "sample_multi":
            expected = len(self.declared_samples)
            if n_cols != expected:
                self._err(
                    "E_COL_001",
                    f"sample/multi mode expects {expected} sample columns, got {n_cols}.",
                    sv_id,
                    line_no,
                    blocking=True,
                )
            return

        if self.mode == "caller_merge":
            callers = _parse_sources(info)
            if n_cols != len(callers):
                self._err(
                    "E_COL_002",
                    f"caller-merge mode has {n_cols} caller column(s), "
                    f"but SOURCES lists {len(callers)} caller(s).",
                    sv_id,
                    line_no,
                    blocking=True,
                )
            return

        if self.mode == "single":
            if n_cols != 1:
                self._err(
                    "E_COL_003",
                    f"single-caller mode expects 1 sample/caller column, got {n_cols}.",
                    sv_id,
                    line_no,
                    blocking=True,
                )
            return

        if self.mode == "mixed":
            # Avoid cascaded false assumptions after the global mixed-mode error.
            return

    def _check_svtype_and_coords(
        self,
        info: dict,
        alt: str,
        pos: str,
        sv_id: str,
        line_no: int,
    ) -> None:
        svtype = info.get("SVTYPE")
        if svtype is None:
            return

        if svtype not in LEGAL_SVTYPES:
            self._err(
                "E_SVTYPE_001",
                f"Illegal SVTYPE '{svtype}'. Allowed: {sorted(LEGAL_SVTYPES)}.",
                sv_id,
                line_no,
                blocking=True,
            )
            return

        end = info.get("END")

        if svtype in {"DEL", "DUP", "INV"}:
            self._check_span_sv_end(svtype, end, pos, sv_id, line_no)

        elif svtype in {"TRA", "BND"}:
            self._check_tra_bnd(svtype, alt, info, sv_id, line_no)

        # INS: SVTYPE legal + CHROM/POS checked.
        # END/SVLEN/CHR2 keys are required by the SVCF INFO contract above,
        # but values may be '.', and no BND bracket is required.

    def _check_span_sv_end(
        self,
        svtype: str,
        end: object,
        pos: str,
        sv_id: str,
        line_no: int,
    ) -> None:
        """DEL/DUP/INV require numeric END >= POS."""
        if end in (None, ".", True):
            self._err(
                "E_END_001",
                f"{svtype} requires a numeric END.",
                sv_id,
                line_no,
                blocking=True,
            )
            return

        if not str(end).isdigit():
            self._err(
                "E_END_002",
                f"{svtype} END must be an integer, got '{end}'.",
                sv_id,
                line_no,
                blocking=True,
            )
            return

        if _is_positive_int(pos) and int(str(end)) < int(pos):
            self._err(
                "E_END_003",
                f"{svtype} END ({end}) < POS ({pos}).",
                sv_id,
                line_no,
                blocking=True,
            )

    def _check_tra_bnd(
        self,
        svtype: str,
        alt: str,
        info: dict,
        sv_id: str,
        line_no: int,
    ) -> None:
        """TRA/BND require CHR2, numeric END, parseable ALT, and agreement."""
        chr2 = info.get("CHR2")
        end = info.get("END")

        if chr2 in (None, "", ".", True):
            self._err(
                "E_TRA_001",
                f"{svtype} requires CHR2.",
                sv_id,
                line_no,
                blocking=True,
            )

        if end in (None, "", ".", True) or not str(end).isdigit():
            self._err(
                "E_TRA_002",
                f"{svtype} requires a numeric END as mate position.",
                sv_id,
                line_no,
                blocking=True,
            )

        parsed = _parse_bnd_alt(alt)
        if parsed is None:
            self._err(
                "E_TRA_003",
                f"{svtype} ALT is not a valid BND bracket form: '{alt}'.",
                sv_id,
                line_no,
                blocking=True,
            )
            return

        mate_chrom, mate_pos = parsed

        if chr2 not in (None, "", ".", True) and mate_chrom != str(chr2):
            self._err(
                "E_TRA_004",
                f"{svtype} ALT mate chrom '{mate_chrom}' != CHR2 '{chr2}'.",
                sv_id,
                line_no,
                blocking=True,
            )

        if str(end).isdigit() and mate_pos != int(str(end)):
            self._err(
                "E_TRA_005",
                f"{svtype} ALT mate pos {mate_pos} != END {end}.",
                sv_id,
                line_no,
                blocking=True,
            )

    def _check_co(
        self,
        sample_cols: list[str],
        sv_id: str,
        line_no: int,
    ) -> None:
        """Check whether at least one sample/caller column has a parseable CO.

        CO is the final FORMAT field under OctopuSV's SVCF contract. Because
        legal sample/caller values may contain ':' inside IDs or ALT strings,
        we must not validate sample/caller values by raw ':' counting.
        """
        found_parseable = False

        for col in sample_cols:
            co = col.rsplit(":", 1)[-1]
            if co and co != "." and "-" in co:
                found_parseable = True
                break

        if found_parseable:
            return

        if self.strict_co:
            self._err(
                "E_CO_001",
                "No parseable CO token in any sample/caller column.",
                sv_id,
                line_no,
                blocking=False,
            )
        else:
            self._warn(
                "W_CO_001",
                "No parseable CO token in any sample/caller column.",
                sv_id,
                line_no,
            )

    # ------------------------------------------------------------------
    # Reporting
    # ------------------------------------------------------------------

    @property
    def errors(self) -> list[Issue]:
        return [issue for issue in self.issues if issue.level == "error"]

    @property
    def warnings(self) -> list[Issue]:
        return [issue for issue in self.issues if issue.level == "warning"]

    @property
    def blocking(self) -> bool:
        return any(issue.blocking for issue in self.issues)

    def status(self) -> str:
        if self.unreadable:
            return "UNREADABLE"
        if self.errors:
            return "FAILED"
        if self.warnings:
            return "PASS_WITH_WARNINGS"
        return "PASS"

    def exit_code(self) -> int:
        if self.unreadable:
            return 2
        if self.errors:
            return 1
        return 0

    def _limited_issues(self, max_issues: int | None) -> tuple[list[Issue], int]:
        if max_issues is None or max_issues <= 0:
            return self.issues, 0

        shown = self.issues[:max_issues]
        truncated = max(0, len(self.issues) - len(shown))
        return shown, truncated

    def to_summary(self, max_issues: int | None = 50) -> str:
        """Human-readable multi-line summary string."""
        lines = [
            "SVCF validation summary",
            f"Input: {self.path}",
            f"Mode: {self.mode}",
            f"Records: {self.records_total}",
            f"Errors: {len(self.errors)}",
            f"Warnings: {len(self.warnings)}",
            f"Status: {self.status()}",
            f"Blocking for Downstream: {self.blocking}",
        ]

        shown_issues, truncated = self._limited_issues(max_issues)

        shown_errors = [issue for issue in shown_issues if issue.level == "error"]
        shown_warnings = [issue for issue in shown_issues if issue.level == "warning"]

        if shown_errors:
            lines.append("")
            lines.append("Errors:")
            for issue in shown_errors:
                loc = f" (line {issue.line_no})" if issue.line_no else ""
                rid = f" [{issue.record_id}]" if issue.record_id else ""
                lines.append(f"  [{issue.code}]{rid}{loc} {issue.message}")

        if shown_warnings:
            lines.append("")
            lines.append("Warnings:")
            for issue in shown_warnings:
                loc = f" (line {issue.line_no})" if issue.line_no else ""
                rid = f" [{issue.record_id}]" if issue.record_id else ""
                lines.append(f"  [{issue.code}]{rid}{loc} {issue.message}")

        if truncated:
            lines.append("")
            lines.append(f"... {truncated} additional issue(s) not shown.")

        return "\n".join(lines)

    def to_json(self, max_issues: int | None = None) -> str:
        """Structured result for agent consumption."""
        shown_issues, truncated = self._limited_issues(max_issues)

        return json.dumps(
            {
                "input": self.path,
                "mode": self.mode,
                "valid": not self.errors and not self.unreadable,
                "status": self.status(),
                "records_total": self.records_total,
                "errors_count": len(self.errors),
                "warnings_count": len(self.warnings),
                "blocking_for_downstream": self.blocking,
                "issues_count": len(self.issues),
                "issues_truncated": truncated,
                "issues": [issue.to_dict() for issue in shown_issues],
            },
            indent=2,
        )
