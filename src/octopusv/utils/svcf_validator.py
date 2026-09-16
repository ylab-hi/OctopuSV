"""SVCF validator: the OctopuSV SVCF contract checker.

This is not a generic VCF lint.  It validates the SVCF structure used by
OctopuSV in three modes:

* single-caller: one evidence column and no non-empty SOURCES field;
* caller-merge: variable-width evidence columns aligned to SOURCES;
* sample/multi: fixed-width sample columns declared by the #CHROM header.

Validation is streaming: records are checked as they are read rather than
loading the whole file into memory.
"""

from __future__ import annotations

import json
import os
import re
from dataclasses import dataclass

from octopusv.utils.svcf_coordinate_parser import parse_svcf_co
from octopusv.utils.svcf_sample_parser import parse_svcf_sample_block
from octopusv.utils.svcf_schema import (
    CALLER_FORMAT,
    SAMPLE_FORMAT,
    MODE_CALLER,
    MODE_MULTI,
    SVCFIdentity,
    parse_identity_from_meta_lines,
    validate_format_for_identity,
    validate_header_columns_for_identity,
    validate_versioned_identity,
)


CALLER_FORMAT_KEYS = CALLER_FORMAT.split(":")
SAMPLE_FORMAT_KEYS = SAMPLE_FORMAT.split(":")

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

# Matching bracket forms only.  The local replacement sequence is deliberately
# not restricted to A/C/G/T/N here; the structural requirement is the bracketed
# remote chrom:pos coordinate.
BND_ALT_RE = re.compile(
    r"^([^:\[\]]*)([\[\]])([^:\[\]]+):(\d+)\2([^:\[\]]*)$"
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
    """Parse INFO into a dict; flags without '=' map to True."""
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


def _duplicate_info_keys(info_str: str) -> list[str]:
    """Return repeated INFO keys in first-repeat order."""
    seen: set[str] = set()
    duplicates: list[str] = []

    for item in info_str.split(";"):
        if not item:
            continue
        key = item.split("=", 1)[0]
        if not key:
            continue
        if key in seen and key not in duplicates:
            duplicates.append(key)
        seen.add(key)

    return duplicates


def _has_sources(info: dict) -> bool:
    value = info.get("SOURCES")
    return value not in (None, "", ".", True)


def _parse_sources(info: dict) -> list[str]:
    """Parse SOURCES in order; duplicate source labels are intentionally kept."""
    if not _has_sources(info):
        return []
    return str(info["SOURCES"]).split(",")


def _source_ids_present(info: dict) -> bool:
    """SOURCE_IDS='.' is present and represents one positional missing ID."""
    return "SOURCE_IDS" in info


def _parse_source_ids(info: dict) -> list[str]:
    """Parse SOURCE_IDS positionally, preserving '.' placeholders."""
    if not _source_ids_present(info):
        return []

    value = info.get("SOURCE_IDS")
    if value is True or value == "":
        return []

    return str(value).split(",")


def _is_positive_int(value: object) -> bool:
    return str(value).isdigit() and int(str(value)) > 0


def _is_nonnegative_int(value: object) -> bool:
    return str(value).isdigit() and int(str(value)) >= 0


def _parse_bnd_alt(alt: str) -> tuple[str, int] | None:
    match = BND_ALT_RE.fullmatch(alt)
    if not match:
        return None

    _, _, mate_chrom, mate_pos, _ = match.groups()
    try:
        return mate_chrom, int(mate_pos)
    except ValueError:
        return None


class SVCFValidator:
    """Validate one SVCF file against the current OctopuSV SVCF contract."""

    def __init__(self, path: str, strict_co: bool = False):
        self.path = path
        self.strict_co = strict_co

        self.issues: list[Issue] = []
        self.records_total = 0
        self.mode: str | None = None
        self.declared_samples: list[str] = []
        self.unreadable = False

        self._has_multi_marker = False
        self.svcf_version: str | None = None
        self.declared_mode: str | None = None
        self._meta_lines: list[str] = []
        self._mode_source_flag: bool | None = None
        self._mixed_mode_reported = False
        self._mode_observations = 0
        self._identity: SVCFIdentity | None = None

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
        """Run all checks in one streaming pass and populate ``issues``."""
        # Make repeat calls deterministic rather than accumulating old state.
        self.issues = []
        self.records_total = 0
        self.mode = None
        self.declared_samples = []
        self.unreadable = False
        self._has_multi_marker = False
        self.svcf_version: str | None = None
        self.declared_mode: str | None = None
        self._meta_lines: list[str] = []
        self._mode_source_flag = None
        self._mixed_mode_reported = False
        self._mode_observations = 0
        self._identity = None

        if not os.path.exists(self.path) or os.path.getsize(self.path) == 0:
            self.unreadable = True
            self._err("E_FILE_001", f"File not found or empty: {self.path}")
            return

        chrom_line: str | None = None
        header_complete = False

        try:
            with open(self.path, encoding="utf-8-sig") as handle:
                for line_no, raw_line in enumerate(handle, start=1):
                    line = raw_line.rstrip("\r\n")
                    if not line:
                        continue

                    if not header_complete:
                        if line.startswith("##"):
                            self._meta_lines.append(line)
                            continue

                        if line.startswith("#CHROM"):
                            chrom_line = line
                            self._initialize_declared_identity()
                            self._check_header(chrom_line)
                            self._initialize_mode_from_header(chrom_line)
                            header_complete = True
                            continue

                        if line.startswith("#"):
                            continue

                        # Data before #CHROM is not a legal SVCF layout.  Stop
                        # instead of trying to infer a header after data began.
                        self._err(
                            "E_HDR_001",
                            "Missing #CHROM header line before data records.",
                            line_no=line_no,
                            blocking=True,
                        )
                        return

                    if line.startswith("#"):
                        continue

                    self.records_total += 1
                    self._check_record(line, line_no)

        except (OSError, UnicodeDecodeError) as exc:
            self.unreadable = True
            self._err("E_FILE_002", f"Cannot read file: {exc}")
            return

        if chrom_line is None:
            self._err(
                "E_HDR_001",
                "Missing #CHROM header line.",
                blocking=True,
            )
            return

        if not self._has_multi_marker and self._mode_observations == 0:
            self.mode = "single"
            if self.records_total == 0:
                self._warn(
                    "W_MODE_001",
                    "No data records found; treating file as single-caller mode.",
                )

    # ------------------------------------------------------------------
    # Header / mode
    # ------------------------------------------------------------------

    def _initialize_declared_identity(self) -> None:
        """Parse explicit SVCF version/mode declarations once per file."""
        try:
            identity = parse_identity_from_meta_lines(self._meta_lines)
        except ValueError as exc:
            self._err(
                "E_HDR_003",
                str(exc),
                blocking=True,
            )
            return

        self._identity = identity
        self.svcf_version = identity.version
        self.declared_mode = identity.mode
        self._has_multi_marker = identity.mode == MODE_MULTI

        try:
            validate_versioned_identity(identity)
        except ValueError as exc:
            self._err(
                "E_VER_001",
                str(exc),
                blocking=True,
            )

    def _check_header(self, chrom_line: str | None) -> None:
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

    def _initialize_mode_from_header(self, chrom_line: str) -> None:
        columns = chrom_line.split("\t")
        trailing_columns = columns[9:]

        if self._identity is not None and self._identity.is_versioned:
            try:
                validate_header_columns_for_identity(
                    self._identity,
                    len(trailing_columns),
                )
            except ValueError as exc:
                self._err(
                    "E_MODE_001",
                    str(exc),
                    blocking=True,
                )

        if self.declared_mode == MODE_MULTI or self._has_multi_marker:
            self.mode = "sample_multi"
            self.declared_samples = trailing_columns
            if not self.declared_samples:
                self._err(
                    "E_MODE_001",
                    "sample/multi mode but #CHROM declares no sample columns.",
                    blocking=True,
                )

    def _update_nonmulti_mode(self, info: dict) -> None:
        """Update no-marker mode without storing prior records."""
        if self.declared_mode == MODE_MULTI or self._has_multi_marker:
            return

        has_sources = _has_sources(info)
        self._mode_observations += 1

        if self._mode_source_flag is None:
            self._mode_source_flag = has_sources
            self.mode = "caller_merge" if has_sources else "single"
            return

        if has_sources == self._mode_source_flag:
            return

        self.mode = "mixed"
        if not self._mixed_mode_reported:
            self._err(
                "E_MODE_002",
                "Invalid mixed SVCF mode: some records have SOURCES and some do not. "
                "A no-marker SVCF must be either all single-caller or all caller-merge.",
                blocking=True,
            )
            self._mixed_mode_reported = True

    # ------------------------------------------------------------------
    # Per-record checks
    # ------------------------------------------------------------------

    def _check_record(self, line: str, line_no: int) -> None:
        parts = line.split("\t")

        # Mode inference needs INFO when available, even if a later column-count
        # check will reject the row.
        info = _parse_info(parts[7]) if len(parts) >= 8 else {}
        if len(parts) >= 8:
            self._update_nonmulti_mode(info)

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

        duplicates = _duplicate_info_keys(info_str)
        if duplicates:
            self._err(
                "E_INFO_002",
                f"Duplicate INFO field(s): {duplicates}.",
                sv_id,
                line_no,
                blocking=True,
            )

        self._check_core_fields(chrom, pos, sv_id, line_no)
        self._check_format(fmt, sv_id, line_no)
        self._check_evidence_width(fmt, sample_cols, sv_id, line_no)
        self._check_info_keys(info, sv_id, line_no)
        self._check_support(info, sv_id, line_no)
        self._check_sources_by_mode(info, sv_id, line_no)
        self._check_column_count(info, sample_cols, sv_id, line_no)
        self._check_source_ids(info, fmt, sample_cols, sv_id, line_no)
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

    def _expected_format(self) -> str:
        """Return the fixed FORMAT schema for the declared/inferred mode."""
        if self.declared_mode == MODE_MULTI or self.mode == "sample_multi":
            return SAMPLE_FORMAT
        return CALLER_FORMAT

    def _expected_format_keys(self) -> list[str]:
        if self.declared_mode == MODE_MULTI or self.mode == "sample_multi":
            return SAMPLE_FORMAT_KEYS
        return CALLER_FORMAT_KEYS

    def _check_format(self, fmt: str, sv_id: str, line_no: int) -> None:
        if self._identity is not None and self._identity.is_versioned:
            try:
                validate_format_for_identity(self._identity, fmt)
            except ValueError as exc:
                self._err(
                    "E_FMT_001",
                    str(exc),
                    sv_id,
                    line_no,
                    blocking=True,
                )
            return

        expected = self._expected_format()
        if fmt != expected:
            self._err(
                "E_FMT_001",
                f"FORMAT for {self.mode or 'unknown'} mode must be '{expected}', "
                f"got '{fmt}'.",
                sv_id,
                line_no,
                blocking=True,
            )

    def _check_evidence_width(
        self,
        fmt: str,
        sample_cols: list[str],
        sv_id: str,
        line_no: int,
    ) -> None:
        """Catch definitely truncated evidence blocks without raw ':' guessing.

        Colon-bearing ID/ALT values may increase the raw token count, but a
        valid fixed SVCF block cannot contain fewer tokens than FORMAT keys.
        """
        expected = self._expected_format()
        if fmt != expected:
            return

        minimum_tokens = len(self._expected_format_keys())
        for index, block in enumerate(sample_cols, start=1):
            token_count = block.count(":") + 1
            if token_count < minimum_tokens:
                self._err(
                    "E_FMT_002",
                    f"Evidence column {index} has only {token_count} colon-delimited "
                    f"token(s); fixed SVCF FORMAT requires at least {minimum_tokens}.",
                    sv_id,
                    line_no,
                    blocking=True,
                )

    def _check_info_keys(self, info: dict, sv_id: str, line_no: int) -> None:
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
        support = info.get("SUPPORT")

        if support is None or support == ".":
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
            sources = _parse_sources(info)
            if n_cols != len(sources):
                self._err(
                    "E_COL_002",
                    f"caller-merge mode has {n_cols} evidence column(s), "
                    f"but SOURCES lists {len(sources)} item(s).",
                    sv_id,
                    line_no,
                    blocking=True,
                )
            return

        if self.mode == "single":
            if n_cols != 1:
                self._err(
                    "E_COL_003",
                    f"single-caller mode expects 1 evidence column, got {n_cols}.",
                    sv_id,
                    line_no,
                    blocking=True,
                )

    def _check_source_ids(
        self,
        info: dict,
        fmt: str,
        sample_cols: list[str],
        sv_id: str,
        line_no: int,
    ) -> None:
        """Check positional SOURCE_IDS without inventing source identity.

        SOURCE_IDS remains optional for legacy SVCF until SVCF 1.1 is finalized.
        When present, its positions are meaningful and '.' placeholders must be
        retained.  Caller-mode evidence IDs are checked directly with the shared
        evidence parser.  We do not infer identity from SC, filenames, or ID
        prefixes.
        """
        if not _source_ids_present(info):
            return

        source_ids = _parse_source_ids(info)
        sources = _parse_sources(info)

        if len(source_ids) != len(sources):
            self._err(
                "E_SRC_003",
                f"SOURCE_IDS lists {len(source_ids)} item(s), but SOURCES lists "
                f"{len(sources)} item(s).",
                sv_id,
                line_no,
                blocking=True,
            )
            return

        if self.mode != "caller_merge":
            return

        if len(sample_cols) != len(source_ids):
            # _check_column_count already reports the evidence/SOURCES mismatch.
            return

        for index, (source_id, sample_col) in enumerate(
            zip(source_ids, sample_cols),
            start=1,
        ):
            parsed = parse_svcf_sample_block(fmt, sample_col)
            evidence_id = str(parsed.get("ID", "."))

            if str(source_id) != evidence_id:
                self._err(
                    "E_SRC_004",
                    f"SOURCE_IDS item {index} ('{source_id}') does not match "
                    f"evidence column {index} ID ('{evidence_id}').",
                    sv_id,
                    line_no,
                    blocking=True,
                )

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

    def _check_span_sv_end(
        self,
        svtype: str,
        end: object,
        pos: str,
        sv_id: str,
        line_no: int,
    ) -> None:
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
        chr2 = info.get("CHR2")
        end = info.get("END")
        svlen = info.get("SVLEN")

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
        else:
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

        if svlen is not None and svlen != ".":
            self._err(
                "E_TRA_006",
                f"{svtype} requires SVLEN='.', got '{svlen}'.",
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
        """Require at least one unambiguously parseable CO value."""
        found_parseable = False

        for col in sample_cols:
            co = col.rsplit(":", 1)[-1]
            if parse_svcf_co(co) is not None:
                found_parseable = True
                break

        if found_parseable:
            return

        if self.strict_co:
            self._err(
                "E_CO_001",
                "No unambiguously parseable CO value in any evidence column.",
                sv_id,
                line_no,
                blocking=False,
            )
        else:
            self._warn(
                "W_CO_001",
                "No unambiguously parseable CO value in any evidence column.",
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
        lines = [
            "SVCF validation summary",
            f"Input: {self.path}",
            f"SVCF Version: {self.svcf_version or 'legacy/unversioned'}",
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

        issue_codes = {issue.code for issue in self.issues}
        if issue_codes & {"E_FMT_002", "E_SRC_004"}:
            lines.append("")
            lines.append("Migration note:")
            lines.append(
                "  If this SVCF was generated by an earlier OctopuSV release, "
                "re-run the original merge with OctopuSV 0.5.0 rather than "
                "manually editing SOURCE_IDS or evidence columns."
            )

        return "\n".join(lines)

    def to_json(self, max_issues: int | None = None) -> str:
        shown_issues, truncated = self._limited_issues(max_issues)

        return json.dumps(
            {
                "input": self.path,
                "svcf_version": self.svcf_version,
                "declared_mode": self.declared_mode,
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
