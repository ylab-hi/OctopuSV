"""SVCF/VCF header reader for `octopusv header`.

Reads ONLY the leading header block (lines starting with '#') and stops at the
first data line. It never inspects SV records, so everything it reports is the
*declared contract* of the file -- what the header says, not what the data
actually contains.

Two outputs:
  - raw dump: every '#' line, verbatim, in order.
  - declared-contract JSON: structured header facts for agents.

Field naming is deliberately honest:
  * mode_hint               -- best the header alone can tell; 'sample_multi'
                               when the marker is present, otherwise
                               'single_or_caller_merge' (single vs caller-merge
                               can only be told by scanning records).
  * trailing_columns        -- columns after FORMAT; could be samples or
                               callers depending on mode, so not called
                               'sample_columns'.
  * declared_info_keys /
    declared_format_keys    -- pulled from ##INFO=<ID=..> / ##FORMAT=<ID=..>
                               definition lines; a declared key may be unused
                               in data, and data could carry undeclared keys.
  * inferred_build          -- genome build guessed from ##contig lengths only;
                               conservative ('unknown' unless key contigs agree).
"""

import json

from octopusv.utils.genome_build import infer_build

CORE_COLUMNS = ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"]
MODE_MULTI_MARKER = "##OctopuSV_mode=multi"


def _extract_id(line):
    """Pull the ID from a ##INFO=<ID=X,..> or ##FORMAT=<ID=X,..> line.

    Returns the ID string or None if the line is not in that shape.
    """
    if "<ID=" not in line:
        return None
    try:
        return line.split("<ID=", 1)[1].split(",", 1)[0].split(">", 1)[0]
    except IndexError:
        return None


class HeaderReader:
    """Read and summarize the header block of an SVCF/VCF file."""

    def __init__(self, path):
        self.path = path
        self.header_lines = []        # every '#' line, verbatim
        self.chrom_line = None        # the '#CHROM' line, if present
        self.unreadable = False

    def read(self):
        """Read header lines; stop at the first non-'#' (data) line."""
        try:
            with open(self.path, encoding="utf-8-sig") as fh:
                for line in fh:
                    if line.startswith("#"):
                        stripped = line.rstrip("\n")
                        self.header_lines.append(stripped)
                        if stripped.startswith("#CHROM"):
                            self.chrom_line = stripped
                            # #CHROM is the last header line in VCF/SVCF; the
                            # next line is data, so we can stop here.
                            break
                    else:
                        # Reached data before any #CHROM (malformed but possible).
                        break
        except (OSError, UnicodeDecodeError):
            self.unreadable = True

    # -- raw dump ------------------------------------------------------------

    def to_raw(self):
        """Return all header lines joined verbatim."""
        return "\n".join(self.header_lines)

    # -- declared contract ---------------------------------------------------

    def _meta_lines(self):
        return [ln for ln in self.header_lines if ln.startswith("##")]

    def _trailing_columns(self):
        if self.chrom_line is None:
            return []
        cols = self.chrom_line.split("\t")
        return cols[9:] if len(cols) > 9 else []

    def _declared_keys(self, prefix):
        keys = []
        for ln in self.header_lines:
            if ln.startswith(prefix):
                key = _extract_id(ln)
                if key:
                    keys.append(key)
        return keys

    def _contig_lengths(self):
        """Map contig id -> length parsed from ##contig=<ID=..,length=..> lines.

        Returns a dict of raw contig id -> int length. Used only for build
        inference; ids are normalized inside infer_build.
        """
        out = {}
        for ln in self._meta_lines():
            if not ln.startswith("##contig="):
                continue
            cid = _extract_id(ln)
            if cid is None:
                continue
            if "length=" not in ln:
                continue
            try:
                raw = ln.split("length=", 1)[1].split(",", 1)[0].split(">", 1)[0]
                length = int(raw)
            except (IndexError, ValueError):
                continue
            out[cid] = length
        return out

    def to_contract(self):
        """Return the declared-contract dict (header facts only)."""
        meta = self._meta_lines()
        has_marker = any(ln.strip() == MODE_MULTI_MARKER for ln in meta)
        chrom_cols = self.chrom_line.split("\t") if self.chrom_line else []

        if has_marker:
            mode_hint = "sample_multi"
            requires_scan = False
        else:
            # Without the marker, the header cannot distinguish single from
            # caller-merge: that depends on whether records carry SOURCES.
            mode_hint = "single_or_caller_merge"
            requires_scan = True

        return {
            "input": self.path,
            "has_chrom_header": self.chrom_line is not None,
            "has_octopusv_mode_marker": has_marker,
            "mode_hint": mode_hint,
            "exact_mode_requires_record_scan": requires_scan,
            "core_columns": chrom_cols[:9] if chrom_cols else [],
            "trailing_columns": self._trailing_columns(),
            "n_meta_lines": len(meta),
            "n_contig_definitions": sum(1 for ln in meta if ln.startswith("##contig=")),
            "declared_info_keys": self._declared_keys("##INFO="),
            "declared_format_keys": self._declared_keys("##FORMAT="),
            "inferred_build": infer_build(self._contig_lengths()),
        }

    def to_json(self):
        return json.dumps(self.to_contract(), indent=2)