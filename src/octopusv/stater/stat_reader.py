"""Single-pass SVCF reader shared by all stat analyzers.

The old stat code opened the file once per analyzer (5x IO on large files) and
each analyzer parsed INFO its own (fragile) way. This module reads the file
ONCE, parsing every record into a small structured dict, and every analyzer
consumes that same in-memory list. INFO parsing is centralized and safe
(handles flags and values containing '=').
"""

from octopusv.utils.svcf_schema import (
    MODE_MULTI,
    parse_identity_from_meta_lines,
    validate_format_for_identity,
    validate_header_columns_for_identity,
    validate_versioned_identity,
)


def parse_info(info_str):
    """Parse an INFO column into a dict. Flags (no '=') map to True.

    Safe against values that themselves contain '=' (splits on the first only),
    unlike dict(item.split('=')).
    """
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


class SVRecord:
    """A minimal structured view of one SVCF data row.

    Only the fields the analyzers need are kept. Sample columns are retained
    raw (one string per caller/sample) so genotype resolution can run later.
    """

    __slots__ = ("chrom", "pos", "sv_id", "qual", "filter", "info",
                 "format", "sample_cols", "svtype")

    def __init__(self, fields):
        self.chrom = fields[0]
        self.pos = fields[1]
        self.sv_id = fields[2]
        self.qual = fields[5]
        self.filter = fields[6]
        self.info = parse_info(fields[7])
        self.format = fields[8] if len(fields) > 8 else ""
        self.sample_cols = fields[9:] if len(fields) > 9 else []
        self.svtype = self.info.get("SVTYPE", "Unknown")


def read_records(input_file):
    """Read an SVCF file once. Returns (records, sample_names, mode).

    Versioned SVCF uses the shared schema contract: unsupported versions and
    mode/FORMAT conflicts fail loudly. Unversioned files retain the historical
    compatibility behavior, including multi-sample inference from either the
    legacy ``##OctopuSV_mode=multi`` marker or multiple trailing columns.
    """
    records = []
    sample_names = []
    meta_lines = []
    identity = None

    with open(input_file) as fh:
        for line in fh:
            stripped = line.rstrip("\r\n")

            if line.startswith("##"):
                meta_lines.append(stripped)
                continue

            if line.startswith("#CHROM"):
                try:
                    identity = parse_identity_from_meta_lines(meta_lines)
                    validate_versioned_identity(identity)
                except ValueError as exc:
                    raise ValueError(
                        f"Invalid SVCF identity in {str(input_file)!r}: {exc}"
                    ) from exc

                header = stripped.split("\t")
                sample_names = header[9:] if len(header) > 9 else []
                try:
                    validate_header_columns_for_identity(
                        identity,
                        len(sample_names),
                    )
                except ValueError as exc:
                    raise ValueError(
                        f"Invalid SVCF header in {str(input_file)!r}: {exc}"
                    ) from exc
                continue

            if line.startswith("#"):
                continue

            fields = stripped.split("\t")
            if len(fields) < 8:
                continue

            if identity is None:
                try:
                    identity = parse_identity_from_meta_lines(meta_lines)
                    validate_versioned_identity(identity)
                except ValueError as exc:
                    raise ValueError(
                        f"Invalid SVCF identity in {str(input_file)!r}: {exc}"
                    ) from exc

            record = SVRecord(fields)
            if identity.is_versioned:
                try:
                    validate_format_for_identity(identity, record.format)
                except ValueError as exc:
                    raise ValueError(
                        f"Invalid SVCF schema in {str(input_file)!r}: {exc}"
                    ) from exc

            records.append(record)

    if identity is None:
        try:
            identity = parse_identity_from_meta_lines(meta_lines)
            validate_versioned_identity(identity)
        except ValueError as exc:
            raise ValueError(
                f"Invalid SVCF identity in {str(input_file)!r}: {exc}"
            ) from exc

    if identity.is_versioned:
        # Once a file declares SVCF 1.1 identity, mode is explicit and column
        # count must not override it. Legacy inference is reserved for
        # unversioned inputs only.
        mode = "sample" if identity.mode == MODE_MULTI else "caller"
    else:
        mode = (
            "sample"
            if identity.mode == MODE_MULTI or len(sample_names) > 1
            else "caller"
        )
    return records, sample_names, mode
