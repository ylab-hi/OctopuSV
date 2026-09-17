"""Shared SVCF schema identity and version contract.

This module deliberately contains only facts that all SVCF producers and
consumers must agree on.  It does not infer caller identity, biological sample
state, or merge topology.
"""

from __future__ import annotations

from dataclasses import dataclass


SVCF_VERSION = "1.1"

MODE_CALLER = "caller"
MODE_MULTI = "multi"
VALID_MODES = {MODE_CALLER, MODE_MULTI}

CALLER_FORMAT = "GT:AD:LN:ST:QV:TY:ID:SC:REF:ALT:CO"
SAMPLE_FORMAT = "GT:AD:UC:UV:LN:ST:QV:TY:ID:SC:REF:ALT:CO"

VERSION_PREFIX = "##SVCFVersion="
MODE_PREFIX = "##OctopuSV_mode="


@dataclass(frozen=True)
class SVCFIdentity:
    """Explicit SVCF identity declared by meta-header lines."""

    version: str | None = None
    mode: str | None = None

    @property
    def is_versioned(self) -> bool:
        return self.version is not None

    @property
    def is_v11(self) -> bool:
        return self.version == SVCF_VERSION


def version_header(version: str = SVCF_VERSION) -> str:
    return f"{VERSION_PREFIX}{version}"


def mode_header(mode: str) -> str:
    if mode not in VALID_MODES:
        raise ValueError(
            f"Invalid SVCF mode {mode!r}; expected one of {sorted(VALID_MODES)}."
        )
    return f"{MODE_PREFIX}{mode}"




def validate_positional_info_item(
    value,
    *,
    field_name: str,
    allow_dot: bool = False,
) -> str:
    """Validate one atom used in SVCF positional INFO lists.

    SVCF 1.1 uses commas to separate positional items inside SOURCES and
    SOURCE_IDS, while semicolons and equals signs delimit INFO fields.  The
    format intentionally defines no escaping layer for these atoms.

    ``.`` is allowed only where the caller explicitly permits a positional
    missing value (SOURCE_IDS).  Known identities must be represented
    explicitly rather than by an ambiguous placeholder.
    """
    text = str(value)

    if text == "":
        raise ValueError(f"{field_name} item must not be empty.")

    if text == ".":
        if allow_dot:
            return text
        raise ValueError(
            f"{field_name} item '.' is reserved for missing values and "
            "cannot be used as a source identity."
        )

    if any(ch in text for ch in ",;=") or any(ch.isspace() for ch in text):
        raise ValueError(
            f"{field_name} item {text!r} cannot be represented safely in "
            "SVCF 1.1 positional INFO lists; items must not contain ',', "
            "';', '=', or whitespace."
        )

    return text


def validate_source_label(value) -> str:
    """Validate one SVCF source label used by SOURCES / FORMAT SC.

    ``:`` is reserved by the fixed evidence-block encoding.  Unlike record IDs
    and ALT, source labels have no escaping layer and therefore cannot contain
    a colon without making ``ID:SC:REF:ALT:CO`` ambiguous.
    """
    text = validate_positional_info_item(
        value,
        field_name="SOURCES",
        allow_dot=False,
    )
    if ":" in text:
        raise ValueError(
            f"SOURCES item {text!r} cannot contain ':' because SVCF 1.1 "
            "uses ':' to delimit FORMAT evidence fields, including SC."
        )
    return text


def validate_source_id(value) -> str:
    """Validate one SOURCE_IDS item; ``.`` is a legal positional placeholder."""
    return validate_positional_info_item(
        value,
        field_name="SOURCE_IDS",
        allow_dot=True,
    )


def parse_identity_from_meta_lines(meta_lines) -> SVCFIdentity:
    """Parse SVCF version/mode declarations from meta-header lines.

    Duplicate declarations with different values are rejected because a file
    cannot have two simultaneous schema identities.
    """

    version = None
    mode = None

    for raw_line in meta_lines:
        line = str(raw_line).strip()

        if line.startswith(VERSION_PREFIX):
            value = line[len(VERSION_PREFIX):].strip()
            if version is not None and value != version:
                raise ValueError(
                    "Conflicting SVCFVersion declarations: "
                    f"{version!r} and {value!r}."
                )
            version = value
            continue

        if line.startswith(MODE_PREFIX):
            value = line[len(MODE_PREFIX):].strip()
            if mode is not None and value != mode:
                raise ValueError(
                    "Conflicting OctopuSV_mode declarations: "
                    f"{mode!r} and {value!r}."
                )
            mode = value

    return SVCFIdentity(version=version, mode=mode)


def expected_format_for_mode(mode: str) -> str:
    if mode == MODE_CALLER:
        return CALLER_FORMAT
    if mode == MODE_MULTI:
        return SAMPLE_FORMAT
    raise ValueError(
        f"Invalid SVCF mode {mode!r}; expected one of {sorted(VALID_MODES)}."
    )


def classify_format(format_field: str) -> str | None:
    """Return the SVCF schema mode represented by one FORMAT string."""

    if format_field == CALLER_FORMAT:
        return MODE_CALLER
    if format_field == SAMPLE_FORMAT:
        return MODE_MULTI
    return None


def validate_versioned_identity(identity: SVCFIdentity) -> None:
    """Validate explicit version/mode declarations.

    Unversioned files are intentionally left to legacy compatibility logic.
    Once a file declares an SVCF version, however, its identity must be
    complete and unambiguous.
    """

    if not identity.is_versioned:
        return

    if identity.version != SVCF_VERSION:
        raise ValueError(
            f"Unsupported SVCFVersion={identity.version!r}; "
            f"this OctopuSV build supports SVCF {SVCF_VERSION}."
        )

    if identity.mode is None:
        raise ValueError(
            f"SVCF {SVCF_VERSION} requires an explicit ##OctopuSV_mode="
            "caller|multi declaration."
        )

    if identity.mode not in VALID_MODES:
        raise ValueError(
            f"Invalid ##OctopuSV_mode={identity.mode!r}; "
            f"expected one of {sorted(VALID_MODES)}."
        )


def validate_header_columns_for_identity(
    identity: SVCFIdentity,
    trailing_column_count: int,
) -> None:
    """Validate the #CHROM trailing-column shape for versioned SVCF.

    SVCF 1.1 caller mode uses one placeholder header column (``SAMPLE``),
    while individual records may carry one or more caller-evidence blocks.
    Multi mode is a biological-sample matrix and therefore requires at least
    one declared sample column.  Unversioned files stay on the legacy path.
    """
    validate_versioned_identity(identity)

    if not identity.is_versioned:
        return

    if identity.mode == MODE_CALLER and trailing_column_count != 1:
        raise ValueError(
            f"SVCF {identity.version} mode 'caller' requires exactly one "
            "#CHROM trailing column; "
            f"got {trailing_column_count}."
        )

    if identity.mode == MODE_MULTI and trailing_column_count < 1:
        raise ValueError(
            f"SVCF {identity.version} mode 'multi' requires at least one "
            "sample column in #CHROM."
        )


def validate_format_for_identity(
    identity: SVCFIdentity,
    format_field: str,
) -> str | None:
    """Validate one record FORMAT against an explicit SVCF identity.

    Unversioned files remain in the legacy compatibility path: their FORMAT is
    classified when possible but is not forced into the 1.1 contract.  Once a
    file declares an SVCF version, however, every record must use the exact
    schema associated with the declared mode.

    Returns the classified schema mode (``caller`` / ``multi``) or ``None``
    for an unversioned, non-SVCF FORMAT.
    """
    validate_versioned_identity(identity)
    schema_mode = classify_format(format_field)

    if not identity.is_versioned:
        return schema_mode

    expected = expected_format_for_mode(identity.mode)
    if format_field != expected:
        raise ValueError(
            f"SVCF {identity.version} mode {identity.mode!r} requires "
            f"FORMAT={expected!r}; got {format_field!r}."
        )

    return schema_mode

