"""SVCF single-record inspection engine.

`octopusv inspect` takes one or more SV IDs and renders each record's PARSED
structure: endpoints, span, provenance (SOURCES/SOURCE_IDS), per-caller or
per-sample evidence blocks, and the raw REF/ALT/INFO.

Boundary (hard):
    inspect explains the SVCF RECORD only. It does NOT annotate genes, predict
    consequences, or make any biological interpretation. Those belong to the
    downstream interpretation layer.

Mode model (hybrid, file_layout + record_mode):
    file_layout  -- structural shape from the header:
        "sample_multi"                 : ##OctopuSV_mode=multi present, or the
                                         #CHROM line has >1 trailing column.
        "single_column_or_caller_merge": one trailing column; cannot tell single
                                         from caller-merge by the header alone.
    record_mode  -- evidence semantics of THIS record:
        "sample" : sample_multi layout -> trailing columns are samples (some may
                   be 0/0 placeholders); identity comes from header sample names.
        "caller" : single-column layout with multiple evidence blocks, OR a
                   single block that carries INFO/SOURCES (a merge product).
                   Identity comes from SOURCES order, zipped with the blocks.
        "single" : single block, no SOURCES (a single-caller corrected record).

Coordinate parsing is delegated to SVCFEvent (svcf_parser); we never re-derive
endpoints here, keeping inspect consistent with the rest of OctopuSV.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Optional

from octopusv.utils.svcf_parser import SVCFFileEventCreator


# SV types with a meaningful linear intra-chromosomal span.
SPAN_SVTYPES = {"DEL", "DUP", "INV"}
NON_LINEAR_SVTYPES = {"TRA", "BND"}

# FORMAT subfields that may themselves contain ':' and therefore must absorb all
# remaining colon-delimited tokens once reached (mirrors SVCFEvent._parse_sample).
_SPECIAL_FORMAT_FIELDS = ("ALT", "CO", "REF")

MODE_MULTI_MARKER = "##OctopuSV_mode=multi"

def parse_format_block(format_keys, block):
    """Split one FORMAT sample/evidence block into a {key: value} dict.

    OctopuSV's standard FORMAT tail is:
        ...:ID:SC:REF:ALT:CO

    ID may contain ':' for some caller/BND-derived IDs.
    ALT may contain ':' for BND/TRA alleles.
    CO is the final field and does not contain ':'.

    Therefore we parse from both ends:
      - fixed fields before ID from the left
      - CO from the far right
      - ALT from the right, allowing BND ALT to contain one ':'
      - REF and SC immediately before ALT
      - everything left in the middle becomes ID
    """
    if not block:
        return {}

    parts = block.split(":")
    result = {}

    # Most OctopuSV SVCF records use this standard tail.
    standard_tail = ["ID", "SC", "REF", "ALT", "CO"]
    has_standard_tail = (
        len(format_keys) >= 5
        and format_keys[-5:] == standard_tail
    )

    if not has_standard_tail:
        # Generic fallback: positional parse; final key absorbs overflow.
        for i, key in enumerate(format_keys):
            if i < len(parts):
                if i == len(format_keys) - 1:
                    result[key] = ":".join(parts[i:])
                else:
                    result[key] = parts[i]
            else:
                result[key] = "."
        return result

    n_head = len(format_keys) - 5

    # Placeholder or malformed short block: positional fill, do not guess.
    if len(parts) < len(format_keys):
        for i, key in enumerate(format_keys):
            result[key] = parts[i] if i < len(parts) else "."
        return result

    # Parse fields before ID from the left.
    for i in range(n_head):
        result[format_keys[i]] = parts[i] if i < len(parts) else "."

    # Remaining tokens correspond to ID:SC:REF:ALT:CO,
    # but ID and ALT may contain ':'.
    tail_parts = parts[n_head:]

    if len(tail_parts) < 5:
        # Should not happen for valid OctopuSV blocks, but keep safe fallback.
        tail_keys = standard_tail
        for i, key in enumerate(tail_keys):
            result[key] = tail_parts[i] if i < len(tail_parts) else "."
        return result

    co_val = tail_parts[-1]
    before_co = tail_parts[:-1]

    # Detect BND/TRA ALT split by ':'.
    # Examples after split:
    #   G[hs37d5:6434738[  -> ["G[hs37d5", "6434738["]
    #   ]10:87115249]T    -> ["]10", "87115249]T"]
    if (
        len(before_co) >= 4
        and ("[" in before_co[-1] or "]" in before_co[-1])
        and ("[" in before_co[-2] or "]" in before_co[-2])
    ):
        alt_tokens = before_co[-2:]
        ref_index = len(before_co) - 3
    else:
        alt_tokens = [before_co[-1]]
        ref_index = len(before_co) - 2

    if ref_index < 1:
        # Not enough tokens for SC + REF + ALT.
        # Fall back conservatively.
        for i, key in enumerate(format_keys):
            result[key] = parts[i] if i < len(parts) else "."
        return result

    ref_val = before_co[ref_index]
    sc_val = before_co[ref_index - 1]
    id_tokens = before_co[:ref_index - 1]
    id_val = ":".join(id_tokens) if id_tokens else "."

    result["ID"] = id_val
    result["SC"] = sc_val
    result["REF"] = ref_val
    result["ALT"] = ":".join(alt_tokens)
    result["CO"] = co_val

    return result


def _split_list_field(value) -> list[str]:
    """Split a comma-separated INFO value into a list; [] for missing/'.'."""
    if value in (None, "", ".", True):
        return []
    return [tok for tok in str(value).split(",") if tok != ""]


def _raw_info_duplicate_keys(raw_info: str) -> list[str]:
    """Return INFO keys appearing more than once in the RAW info string.

    SVCFEvent._parse_info collapses duplicates, so duplicates (e.g. a buggy merge
    writing SOURCE_IDS twice) are invisible in event.info. We scan the raw text.
    """
    if not raw_info:
        return []
    counts: dict[str, int] = {}
    for item in raw_info.split(";"):
        if not item:
            continue
        key = item.split("=", 1)[0]
        counts[key] = counts.get(key, 0) + 1
    return sorted(k for k, n in counts.items() if n > 1)


def _is_placeholder_block(block: str) -> bool:
    """A sample column is a non-supporting placeholder when GT is 0/0 or ./.."""
    if not block:
        return True
    gt = block.split(":", 1)[0]
    return gt in ("0/0", "./.")


class RecordInspector:
    """Build the structured inspection dict for a single SVCFEvent."""

    def __init__(self, event, file_layout: str):
        self.event = event
        self.file_layout = file_layout

    def to_dict(self) -> dict:
        ev = self.event
        warnings: list[dict] = []

        svtype = ev.sv_type or "Unknown"
        format_keys = ev.format.split(":") if ev.format else []

        sources = _split_list_field(ev.info.get("SOURCES"))
        source_ids = _split_list_field(ev.info.get("SOURCE_IDS"))

        # Duplicate INFO keys (from raw INFO, since event.info collapses them).
        for key in _raw_info_duplicate_keys(getattr(ev, "raw_info", "")):
            warnings.append({
                "code": "W_DUPLICATE_INFO_KEY",
                "message": f"INFO key '{key}' appears more than once in this record.",
            })

        record_mode = self._infer_record_mode(sources)

        # SOURCES vs SOURCE_IDS length (provenance pairing sanity).
        if sources and source_ids and len(sources) != len(source_ids):
            warnings.append({
                "code": "W_SOURCE_IDS_LENGTH_MISMATCH",
                "message": (
                    f"SOURCES (n={len(sources)}) and SOURCE_IDS "
                    f"(n={len(source_ids)}) have different lengths."
                ),
            })

        # --- geometry (delegated to SVCFEvent) ------------------------------
        is_interchromosomal = ev.start_chrom != ev.end_chrom
        endpoint1 = {"contig": ev.start_chrom, "pos": ev.start_pos}
        endpoint2 = {"contig": ev.end_chrom, "pos": ev.end_pos}
        span, span_reason = self._compute_span(svtype, ev, is_interchromosomal)
        contigs_involved = sorted({ev.start_chrom, ev.end_chrom})

        # --- evidence -------------------------------------------------------
        if record_mode == "sample":
            evidence, carrier_samples = self._sample_evidence(format_keys)
            # carrier columns should equal INFO/SOURCES.
            if sources and sorted(carrier_samples) != sorted(sources):
                warnings.append({
                    "code": "W_SAMPLE_SOURCES_CARRIER_MISMATCH",
                    "message": (
                        "Carrier sample columns do not match INFO/SOURCES "
                        f"(carriers={carrier_samples}, sources={sources})."
                    ),
                })
            evidence_extra = {"carrier_samples": carrier_samples}
        else:
            evidence, ev_warnings = self._caller_evidence(
                format_keys, sources, source_ids
            )
            warnings.extend(ev_warnings)
            evidence_extra = {}

        record = {
            "id": ev.sv_id,
            "svtype": svtype,
            "file_layout": self.file_layout,
            "record_mode": record_mode,
            "chrom": ev.chrom,
            "pos": ev.pos,
            "chrom2": ev.info.get("CHR2", ev.end_chrom),
            "end": ev.end_pos,
            "endpoint1": endpoint1,
            "endpoint2": endpoint2,
            "span": span,
            "span_reason": span_reason,
            "length_method": "end_minus_start",
            "contigs_involved": contigs_involved,
            "is_interchromosomal": is_interchromosomal,
            "filter": ev.filter,
            "qual": ev.quality,
            "support": ev.info.get("SUPPORT", "."),  # raw, never recomputed
            "sources": sources,
            "source_ids": source_ids,
            "format_keys": format_keys,
            "evidence": evidence,
            "raw": {
                "ref": ev.ref,
                "alt": ev.alt,
                "info": dict(ev.info),  # collapsed; duplicates flagged in warnings
            },
            "warnings": warnings,
        }
        record.update(evidence_extra)
        return record

    def _infer_record_mode(self, sources) -> str:
        """Decide this record's evidence semantics (see module docstring)."""
        if self.file_layout == "sample_multi":
            return "sample"
        n_blocks = len(self.event.raw_sample_columns or [])
        if n_blocks > 1:
            return "caller"
        if sources:
            return "caller"
        return "single"

    def _compute_span(self, svtype, ev, is_interchromosomal):
        if is_interchromosomal or svtype in NON_LINEAR_SVTYPES:
            return None, "interchromosomal_no_linear_span"
        if svtype == "INS":
            return None, "INS_is_endpoint_event_no_reference_span"
        if svtype not in SPAN_SVTYPES:
            return None, f"span_not_defined_for_{svtype}"
        start = min(ev.start_pos, ev.end_pos)
        end = max(ev.start_pos, ev.end_pos)
        return (
            {"contig": ev.start_chrom, "start": start, "end": end, "length": end - start},
            None,
        )

    def _caller_evidence(self, format_keys, sources, source_ids):
        """Caller-mode: zip SOURCES / SOURCE_IDS / evidence blocks (all equal len).

        Identity comes from SOURCES order, which OctopuSV merge guarantees to
        match the evidence-block order.
        """
        blocks = self.event.raw_sample_columns or []
        warnings = []

        lengths = {len(sources), len(source_ids), len(blocks)}
        # Treat empty SOURCE_IDS as "not provided" rather than a length error.
        if source_ids == []:
            lengths = {len(sources), len(blocks)}
        if len(lengths) > 1:
            warnings.append({
                "code": "W_CALLER_EVIDENCE_LENGTH_MISMATCH",
                "message": (
                    "SOURCES, SOURCE_IDS, and caller evidence blocks have "
                    f"different lengths (sources={len(sources)}, "
                    f"source_ids={len(source_ids)}, blocks={len(blocks)})."
                ),
            })

        evidence = []
        n = max(len(sources), len(blocks))
        for i in range(n):
            source = sources[i] if i < len(sources) else None
            sid = source_ids[i] if i < len(source_ids) else None
            block = blocks[i] if i < len(blocks) else None
            fields = parse_format_block(format_keys, block) if block is not None else {}
            evidence.append({
                "index": i,
                "source": source,
                "source_id": sid,
                "is_placeholder": block is not None and _is_placeholder_block(block),
                "fields": fields,
            })
        return evidence, warnings

    def _sample_evidence(self, format_keys):
        """Sample-mode: align each trailing column to its header sample name.

        Identity comes from header sample_names (NOT from SOURCES). Placeholder
        columns (GT 0/0) are kept but flagged; carriers are listed for a
        consistency check against INFO/SOURCES.
        """
        ev = self.event
        sample_names = ev.sample_names or []
        blocks = ev.raw_sample_columns or []

        evidence = []
        carrier_samples = []
        n = max(len(sample_names), len(blocks))
        for i in range(n):
            name = sample_names[i] if i < len(sample_names) else f"col{i}"
            block = blocks[i] if i < len(blocks) else None
            placeholder = block is None or _is_placeholder_block(block)
            fields = parse_format_block(format_keys, block) if block is not None else {}
            if not placeholder:
                carrier_samples.append(name)
            evidence.append({
                "index": i,
                "sample": name,
                "is_carrier": not placeholder,
                "is_placeholder": placeholder,
                "fields": fields,
            })
        return evidence, carrier_samples


class SVCFInspector:
    """Locate records by ID and inspect them, reusing the canonical parser."""

    def __init__(self, input_file: str, ids: list[str]):
        self.input_file = input_file
        self.requested_ids = list(ids)

    def _detect_file_layout(self) -> str:
        """Header-only layout detection."""
        with open(self.input_file, encoding="utf-8-sig") as fh:
            for line in fh:
                if line.startswith("##"):
                    if line.strip() == MODE_MULTI_MARKER:
                        return "sample_multi"
                    continue
                if line.startswith("#CHROM"):
                    cols = line.rstrip("\n").split("\t")
                    if len(cols) > 10:  # >1 trailing column after FORMAT
                        return "sample_multi"
                    return "single_column_or_caller_merge"
                break
        return "single_column_or_caller_merge"

    def run(self) -> dict:
        """Inspect all requested IDs; return a structured result dict.

        NOTE: parses the whole file via SVCFFileEventCreator (the canonical
        parser), then selects by ID. Reads more than strictly needed for huge
        files; kept simple and correct for v0.1. Optimize to an early-exit ID
        scan only if profiling shows it matters.
        """
        if not Path(self.input_file).exists():
            raise FileNotFoundError(f"Input file not found: {self.input_file}")

        file_layout = self._detect_file_layout()

        creator = SVCFFileEventCreator([self.input_file])
        creator.parse()

        by_id: dict[str, object] = {}
        duplicate_ids: set[str] = set()
        for ev in creator.events:
            if ev.sv_id in by_id:
                duplicate_ids.add(ev.sv_id)
                continue
            by_id[ev.sv_id] = ev

        results = []
        for qid in self.requested_ids:
            ev = by_id.get(qid)
            if ev is None:
                results.append({"query_id": qid, "found": False, "record": None})
                continue
            record = RecordInspector(ev, file_layout).to_dict()
            if qid in duplicate_ids:
                record["warnings"].append({
                    "code": "W_DUPLICATE_RECORD_ID",
                    "message": f"ID '{qid}' appears multiple times; showing the first.",
                })
            results.append({"query_id": qid, "found": True, "record": record})

        return {
            "tool": "octopusv inspect",
            "input": self.input_file,
            "file_layout": file_layout,
            "requested": len(self.requested_ids),
            "found": sum(1 for r in results if r["found"]),
            "results": results,
        }

    # -- rendering -----------------------------------------------------------

    @staticmethod
    def to_json(result: dict) -> str:
        return json.dumps(result, indent=2)

    @staticmethod
    def to_text(result: dict) -> str:
        lines = []
        for r in result["results"]:
            qid = r["query_id"]
            if not r["found"]:
                lines.append(f"Record {qid}: NOT FOUND")
                lines.append("")
                continue

            rec = r["record"]
            lines.append("Record")
            lines.append(f"  ID:              {rec['id']}")
            lines.append(f"  SVTYPE:          {rec['svtype']}")
            lines.append(f"  FILTER:          {rec['filter']}")
            lines.append(f"  QUAL:            {rec['qual']}")
            lines.append(f"  file_layout:     {rec['file_layout']}")
            lines.append(f"  record_mode:     {rec['record_mode']}")
            lines.append("")
            lines.append("Coordinates")
            lines.append(f"  endpoint1:       {rec['endpoint1']['contig']}:{rec['endpoint1']['pos']}")
            lines.append(f"  endpoint2:       {rec['endpoint2']['contig']}:{rec['endpoint2']['pos']}")
            if rec["span"] is not None:
                s = rec["span"]
                lines.append(f"  span:            {s['contig']}:{s['start']}-{s['end']}")
                lines.append(f"  length:          {s['length']}")
            else:
                lines.append(f"  span:            n/a ({rec['span_reason']})")
            lines.append(f"  contigs:         {', '.join(rec['contigs_involved'])}")
            lines.append(f"  interchromosomal: {'yes' if rec['is_interchromosomal'] else 'no'}")
            lines.append("")
            lines.append("Provenance")
            lines.append(f"  SOURCES:         {', '.join(rec['sources']) or '.'}")
            lines.append(f"  SOURCE_IDS:      {', '.join(rec['source_ids']) or '.'}")
            lines.append(f"  SUPPORT:         {rec['support']}")
            lines.append("")
            lines.append("Evidence")
            for e in rec["evidence"]:
                label = e.get("source") if "source" in e else e.get("sample")
                tag = " (placeholder)" if e.get("is_placeholder") else ""
                co = e["fields"].get("CO", ".")
                sid = e.get("source_id", "")
                sid_str = f" source_id={sid}" if sid else ""
                lines.append(f"  [{e['index']}] {label}{sid_str} CO={co}{tag}")
            if rec["warnings"]:
                lines.append("")
                lines.append("Warnings")
                for w in rec["warnings"]:
                    lines.append(f"  {w['code']}: {w['message']}")
            lines.append("")

        return "\n".join(lines).rstrip("\n")