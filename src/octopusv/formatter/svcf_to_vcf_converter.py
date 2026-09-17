import io
import logging
import re
from pathlib import Path

from octopusv.utils.atomic_write import atomic_output_path
from octopusv.utils.caller_consensus import resolve_caller_svcf_consensus
from octopusv.utils.svcf_parser import SVCFEvent
from octopusv.utils.svcf_sample_parser import parse_svcf_sample_block
from octopusv.utils.text_io import open_text_auto
from octopusv.utils.svcf_utils import (
    global_meta_audit_lines,
    merge_safe_global_meta_lines,
)
from octopusv.utils.vcf_info import format_vcf_info_item
from octopusv.utils.svcf_schema import (
    MODE_MULTI,
    parse_identity_from_meta_lines,
    validate_format_for_identity,
    validate_header_columns_for_identity,
    validate_versioned_identity,
)
from octopusv.utils.sample_mode_semantics import (
    downstream_sample_gt,
    normalize_unobserved_sample_gt,
)


class SVCFtoVCFConverter:
    """Convert OctopuSV SVCF records back to VCF4.2-compatible records.

    Design goals:
      - Preserve standard VCF compatibility.
      - Preserve OctopuSV source/evidence metadata in INFO, especially
        SOURCES and SOURCE_IDS.
      - Keep sample-mode sample columns intact and in header order.
      - Collapse multi-evidence caller-mode records into one standard VCF
        sample column using the shared order-independent consensus rule.
      - Support streaming conversion so memory use scales with one record,
        not with the complete input/output files.

    ``events`` remains supported for backward compatibility with existing
    programmatic callers and tests.  The CLI uses ``events=None`` and streams
    directly from ``input_svcf_file``.
    """

    DEFAULT_INFO_IDS = {
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
        "SOURCES",
        "SOURCE_IDS",
    }

    DEFAULT_FORMAT_IDS = {
        "GT",
        "AD",
        "DP",
        "UC",
        "UV",
        "LN",
        "ST",
        "QV",
        "TY",
        "ID",
        "SC",
        "REF",
        "ALT",
        "CO",
    }

    STANDARD_INFO_LINES = [
        '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">',
        '##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant described in this record">',
        '##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Difference in length between REF and ALT alleles">',
        '##INFO=<ID=CHR2,Number=1,Type=String,Description="Chromosome for end coordinate">',
        '##INFO=<ID=SUPPORT,Number=1,Type=Integer,Description="Number of reads supporting this variant">',
        '##INFO=<ID=SVMETHOD,Number=1,Type=String,Description="Method used to detect or merge the SV">',
        '##INFO=<ID=RTID,Number=1,Type=String,Description="Associated reciprocal translocation ID if available">',
        '##INFO=<ID=AF,Number=1,Type=Float,Description="Allele frequency if available">',
        '##INFO=<ID=STRAND,Number=1,Type=String,Description="Strand orientation of the SV">',
        '##INFO=<ID=RNAMES,Number=.,Type=String,Description="Supporting read names if available">',
        '##INFO=<ID=SOURCES,Number=.,Type=String,Description="Source callers or samples supporting this converted SVCF record">',
        '##INFO=<ID=SOURCE_IDS,Number=.,Type=String,Description="Original IDs of merged SVs from supporting callers or samples">',
    ]

    STANDARD_ALT_LINES = [
        '##ALT=<ID=DEL,Description="Deletion">',
        '##ALT=<ID=INV,Description="Inversion">',
        '##ALT=<ID=DUP,Description="Duplication">',
        '##ALT=<ID=INS,Description="Insertion">',
        '##ALT=<ID=TRA,Description="Translocation">',
        '##ALT=<ID=BND,Description="Breakend">',
    ]

    STANDARD_FORMAT_LINES = [
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
        '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">',
        '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth derived from AD">',
        '##FORMAT=<ID=LN,Number=1,Type=Integer,Description="Length of SV">',
    ]

    SAMPLE_MODE_FORMAT_LINES = [
        '##FORMAT=<ID=UC,Number=1,Type=Integer,Description="Number of unique callers supporting carrier presence for this synthesized sample call">',
        '##FORMAT=<ID=UV,Number=1,Type=Integer,Description="Number of unique callers contributing a valid presence vote for this synthesized sample call">',
    ]

    def __init__(
        self,
        events=None,
        input_svcf_file=None,
        *,
        unobserved_sample_gt="missing",
    ):
        if input_svcf_file is None:
            raise ValueError("input_svcf_file is required")

        self.events = events
        self.input_svcf_file = str(input_svcf_file)
        self.unobserved_sample_gt = normalize_unobserved_sample_gt(
            unobserved_sample_gt
        )

        (
            self._meta_lines,
            self._contig_lines,
            self._original_definitions,
            self.sample_names,
            self._header_trailing_count,
        ) = self._read_header_metadata()

        try:
            self._identity = parse_identity_from_meta_lines(self._meta_lines)
            validate_versioned_identity(self._identity)
        except ValueError as exc:
            raise ValueError(
                f"Invalid SVCF identity in {self.input_svcf_file!r}: {exc}"
            ) from exc

        try:
            validate_header_columns_for_identity(
                self._identity,
                self._header_trailing_count,
            )
        except ValueError as exc:
            raise ValueError(
                f"Invalid SVCF header in {self.input_svcf_file!r}: {exc}"
            ) from exc

        self._has_multi_marker = self._identity.mode == MODE_MULTI

        if self._identity.is_versioned:
            # Versioned SVCF identity is explicit. Never fall back to column-
            # count inference after a file has declared its schema.
            self.mode = "multi" if self._identity.mode == MODE_MULTI else "single"
        elif self._has_multi_marker:
            # Historical multi-sample files may carry the old mode marker
            # without an SVCFVersion declaration.
            self.mode = "multi"
        elif len(self.sample_names) > 1:
            self.mode = "legacy_multi"
            logging.warning(
                "SVCF has multiple sample columns but no "
                "##OctopuSV_mode=multi marker; treating it as legacy "
                "sample-mode input."
            )
        else:
            self.mode = "single"

        self._first_record_format = self._read_first_record_format()
        if self._first_record_format is not None:
            self._validate_versioned_record_format(self._first_record_format)
        self._supports_unobserved_sample_policy = (
            self.mode == "multi"
            and self._format_has_consensus_fields(self._first_record_format)
        )

        if (
            self.unobserved_sample_gt == "ref"
            and not self._supports_unobserved_sample_policy
        ):
            logging.warning(
                "--unobserved-sample-gt=ref has no effect because this "
                "input does not contain SVCF 1.1 sample consensus fields "
                "(UC/UV)."
            )

    def _read_first_record_format(self) -> str | None:
        """Return FORMAT from the first data record, if one exists.

        The unobserved-sample export policy is meaningful only for SVCF 1.1
        sample-mode records that actually carry UC/UV.  Mode alone is not
        sufficient because legacy multi-sample files can carry the multi
        marker while still using the older 11-field sample schema.
        """
        with open_text_auto(self.input_svcf_file) as handle:
            for raw_line in handle:
                if not raw_line or raw_line.startswith("#"):
                    continue
                parts = raw_line.rstrip("\r\n").split("\t")
                if len(parts) > 8:
                    return parts[8]
                return None
        return None

    @staticmethod
    def _format_has_consensus_fields(format_field: str | None) -> bool:
        if not format_field:
            return False
        keys = set(format_field.split(":"))
        return {"UC", "UV"}.issubset(keys)

    def _validate_versioned_record_format(self, format_field: str) -> None:
        """Enforce the declared SVCF 1.1 schema for one record.

        Unversioned inputs remain on the legacy compatibility path.
        """
        if not self._identity.is_versioned:
            return
        try:
            validate_format_for_identity(self._identity, format_field)
        except ValueError as exc:
            raise ValueError(
                f"Invalid SVCF schema in {self.input_svcf_file!r}: {exc}"
            ) from exc

    @staticmethod
    def _header_id(line: str, kind: str) -> str | None:
        """Extract the ID from a VCF-style header definition line."""
        match = re.search(rf"^##{kind}=<ID=([^,>]+)", line)
        if match:
            return match.group(1)
        return None

    @staticmethod
    def _is_missing_info_value(value) -> bool:
        """Return True for values that should not be emitted in INFO."""
        return value is None or value is True or value == "" or value == "."

    @staticmethod
    def _format_info_value(value) -> str:
        """Make an INFO value safe for VCF text output."""
        text = str(value)
        text = re.sub(r"\s+", "_", text)
        text = text.replace(";", ",")
        text = text.replace("=", "_")
        return text

    def _read_header_metadata(self):
        """Read SVCF header metadata once.

        Returns meta lines, contig lines, preserved original definitions,
        trailing #CHROM sample/evidence-column names, and the exact number of
        trailing header columns.  The explicit count lets versioned SVCF
        validate caller/multi header shape without legacy inference.
        """
        meta_lines = []
        contig_lines = []
        original_definitions = {
            "filter_lines": [],
            "info_lines": [],
            "alt_lines": [],
            "format_lines": [],
        }
        sample_names = ["Sample"]
        header_trailing_count = 0

        with open_text_auto(self.input_svcf_file) as handle:
            for raw_line in handle:
                line = raw_line.rstrip("\r\n")

                if line.startswith("#CHROM"):
                    parts = line.split("\t")
                    header_trailing_count = max(0, len(parts) - 9)
                    if len(parts) > 9:
                        sample_names = parts[9:]
                    break

                if not line.startswith("##"):
                    continue

                meta_lines.append(line)

                if line.startswith("##contig"):
                    contig_lines.append(line)

                elif line.startswith("##FILTER="):
                    filter_id = self._header_id(line, "FILTER")
                    if filter_id != "PASS":
                        original_definitions["filter_lines"].append(line)

                elif line.startswith("##INFO="):
                    info_id = self._header_id(line, "INFO")
                    if info_id is None or info_id not in self.DEFAULT_INFO_IDS:
                        original_definitions["info_lines"].append(line)

                elif line.startswith("##ALT="):
                    alt_id = self._header_id(line, "ALT")
                    if alt_id is None or alt_id not in {
                        "DEL", "DUP", "INV", "INS", "TRA", "BND"
                    }:
                        original_definitions["alt_lines"].append(line)

                elif line.startswith("##FORMAT="):
                    format_id = self._header_id(line, "FORMAT")
                    if format_id is None or format_id not in self.DEFAULT_FORMAT_IDS:
                        original_definitions["format_lines"].append(line)

        return (
            meta_lines,
            contig_lines,
            original_definitions,
            sample_names,
            header_trailing_count,
        )

    def _detect_mode_and_samples(self):
        """Backward-compatible helper returning cached header information."""
        return self.mode, list(self.sample_names)

    def convert(self):
        """Return the complete converted VCF as a string.

        This compatibility API intentionally materializes output in memory.
        The CLI uses ``convert_to_file`` for streaming conversion.
        """
        buffer = io.StringIO()
        self.write(buffer)
        return buffer.getvalue()

    def write(self, handle):
        """Write converted VCF content to an already-open text handle."""
        handle.write(self._generate_vcf_header())
        for event in self._event_iterator():
            handle.write(self._convert_event_to_vcf(event))

    def convert_to_file(self, output_file):
        """Stream conversion to a temporary file, then atomically replace output."""
        with atomic_output_path(output_file) as temp_path:
            with open(temp_path, mode="w", encoding="utf-8") as handle:
                self.write(handle)


    def _event_iterator(self):
        """Return the configured event stream.

        Existing programmatic callers may provide a pre-parsed iterable via
        ``events``.  CLI conversion leaves ``events`` as None and reads one SVCF
        record at a time from disk.
        """
        if self.events is not None:
            return iter(self.events)
        return self._iter_events_from_file()

    def _iter_events_from_file(self):
        """Yield one SVCFEvent at a time without storing the full input file."""
        sample_names = list(self.sample_names)
        sample_name = sample_names[0] if sample_names else "Sample"

        with open_text_auto(self.input_svcf_file) as handle:
            for line_number, line in enumerate(handle, 1):
                if line.startswith("#") or not line.strip():
                    continue

                parts = line.rstrip("\r\n").split("\t")
                if len(parts) < 10:
                    raise ValueError(
                        f"Malformed SVCF record in {str(self.input_svcf_file)!r} "
                        f"on line {line_number}: expected at least 10 "
                        f"tab-separated columns, got {len(parts)}."
                    )

                yield SVCFEvent(
                    *parts[:10],
                    source_file=self.input_svcf_file,
                    sample_name=sample_name,
                    raw_sample_columns=parts[9:],
                    sample_names=sample_names,
                )

    def _extract_original_definitions_from_svcf(self):
        """Return cached non-default header definitions."""
        return {
            key: list(value)
            for key, value in self._original_definitions.items()
        }

    def _read_contig_lines(self) -> list[str]:
        """Return cached ##contig lines from the input SVCF header."""
        return list(self._contig_lines)

    def _append_unique_header_lines(
        self,
        header: str,
        lines: list[str],
        kind: str,
        existing_ids: set[str],
    ) -> str:
        """Append header definition lines while avoiding duplicate IDs."""
        for line in lines:
            line_id = self._header_id(line, kind)
            if line_id is not None and line_id in existing_ids:
                continue
            header += line + "\n"
            if line_id is not None:
                existing_ids.add(line_id)
        return header

    def _generate_vcf_header(self):
        """Generate a VCF4.2 header for the converted VCF."""
        original_defs = self._extract_original_definitions_from_svcf()
        contig_lines = self._read_contig_lines()

        header = "##fileformat=VCFv4.2\n"

        for meta_line in merge_safe_global_meta_lines(
            [(self.input_svcf_file, self._meta_lines)]
        ):
            header += meta_line + "\n"

        # Preserve input-specific reference/assembly audit metadata created by
        # merge when a single standard declaration could not be verified.
        for meta_line in global_meta_audit_lines(self._meta_lines):
            header += meta_line + "\n"

        if self._supports_unobserved_sample_policy:
            header += (
                "##OctopuSV_unobserved_sample_gt="
                f"{self.unobserved_sample_gt}\n"
            )

        for contig_line in contig_lines:
            header += contig_line + "\n"

        existing_filter_ids = {
            self._header_id(line, "FILTER")
            for line in original_defs["filter_lines"]
        }
        existing_filter_ids.discard(None)

        for filter_line in original_defs["filter_lines"]:
            header += filter_line + "\n"

        if "PASS" not in existing_filter_ids:
            header += '##FILTER=<ID=PASS,Description="All filters passed">\n'
            existing_filter_ids.add("PASS")

        existing_info_ids = {
            self._header_id(line, "INFO")
            for line in original_defs["info_lines"]
        }
        existing_info_ids.discard(None)

        for info_line in original_defs["info_lines"]:
            header += info_line + "\n"

        header = self._append_unique_header_lines(
            header,
            self.STANDARD_INFO_LINES,
            "INFO",
            existing_info_ids,
        )

        existing_alt_ids = {
            self._header_id(line, "ALT")
            for line in original_defs["alt_lines"]
        }
        existing_alt_ids.discard(None)

        for alt_line in original_defs["alt_lines"]:
            header += alt_line + "\n"

        header = self._append_unique_header_lines(
            header,
            self.STANDARD_ALT_LINES,
            "ALT",
            existing_alt_ids,
        )

        existing_format_ids = {
            self._header_id(line, "FORMAT")
            for line in original_defs["format_lines"]
        }
        existing_format_ids.discard(None)

        for format_line in original_defs["format_lines"]:
            header += format_line + "\n"

        header = self._append_unique_header_lines(
            header,
            self.STANDARD_FORMAT_LINES,
            "FORMAT",
            existing_format_ids,
        )

        # UC/UV are declared for every converted VCF because caller-mode
        # records with multiple evidence blocks are synthesized into one
        # sample-level call and expose these consensus counts.  Single-evidence
        # caller records keep their original compact FORMAT and simply do not
        # use UC/UV.
        header = self._append_unique_header_lines(
            header,
            self.SAMPLE_MODE_FORMAT_LINES,
            "FORMAT",
            existing_format_ids,
        )

        sample_header = "\t".join(self.sample_names)
        header += (
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"
            f"{sample_header}\n"
        )

        return header

    def _add_info_field(self, info_fields: list[str], event, key: str):
        """Append one INFO field from event.info if it has a meaningful value."""
        value = event.info.get(key)
        if self._is_missing_info_value(value):
            return
        info_fields.append(f"{key}={self._format_info_value(value)}")

    def _convert_event_to_vcf(self, event):
        self._validate_versioned_record_format(event.format)

        chrom = event.chrom
        pos = event.pos
        record_id = event.sv_id

        ref = event.ref if event.ref and event.ref != "." else "N"

        alt = self._get_alt(event)
        qual = event.quality if hasattr(event, "quality") else "."
        filter_status = event.filter if hasattr(event, "filter") else "PASS"

        info_fields = [f"SVTYPE={event.sv_type}"]

        if event.sv_type == "TRA":
            # Bracket-form TRA carries the remote breakpoint in ALT. Symbolic
            # <TRA> does not, so CHR2 + END must remain in the exported VCF or
            # the second breakpoint would be silently lost.
            if alt == "<TRA>":
                end_value = event.info.get("END")
                if not self._is_missing_info_value(end_value):
                    info_fields.append(
                        f"END={self._format_info_value(end_value)}"
                    )
        elif event.sv_type == "BND":
            pass
        elif event.sv_type == "INS":
            info_fields.append(f"END={event.pos}")
        else:
            if hasattr(event, "end_pos") and event.end_pos is not None:
                info_fields.append(f"END={event.end_pos}")

        svlen_value = event.info.get("SVLEN")
        if not self._is_missing_info_value(svlen_value):
            try:
                svlen = int(str(svlen_value).strip())
                if event.sv_type == "DEL":
                    svlen = -abs(svlen)
                elif event.sv_type == "INS":
                    svlen = abs(svlen)

                if event.sv_type != "BND":
                    info_fields.append(f"SVLEN={svlen}")
            except ValueError:
                pass

        for key in [
            "CHR2",
            "SUPPORT",
            "SVMETHOD",
            "RTID",
            "AF",
            "STRAND",
            "RNAMES",
            "SOURCES",
            "SOURCE_IDS",
        ]:
            self._add_info_field(info_fields, event, key)

        # Core conversion-managed fields above are authoritative. Any other
        # event-level INFO annotation that is still representable in ordinary
        # VCF is carried through instead of disappearing silently.
        for key, value in event.info.items():
            if key in self.DEFAULT_INFO_IDS:
                continue
            if value is True:
                info_fields.append(format_vcf_info_item(key, value))
                continue
            if self._is_missing_info_value(value):
                continue
            formatted_value = self._format_info_value(value)
            info_fields.append(format_vcf_info_item(key, formatted_value))

        info = ";".join(info_fields) if info_fields else "."

        raw_cols = getattr(event, "raw_sample_columns", None)
        if not raw_cols:
            format_keys = (event.format or "GT:AD:LN").split(":")
            raw_cols = [
                ":".join(
                    str(event.sample.get(key, "."))
                    for key in format_keys
                )
            ]

        expected = len(self.sample_names)
        is_sample_mode = self.mode in {"multi", "legacy_multi"}

        if is_sample_mode:
            vcf_format = "GT:AD:DP:UC:UV:LN"
            sample_columns = [
                self._convert_svcf_sample_to_vcf(
                    sample,
                    event.format,
                    include_consensus_fields=True,
                )
                for sample in raw_cols
            ]
            if len(sample_columns) != expected:
                raise ValueError(
                    f"Sample column count mismatch for {record_id} at {chrom}:{pos}: "
                    f"got {len(sample_columns)}, expected {expected} "
                    "(from #CHROM header)."
                )
        else:
            # A single caller evidence block is a passthrough, not a synthesis,
            # so preserve the compact historical FORMAT.  Once multiple
            # evidence blocks are collapsed, expose UC/UV so the synthesized
            # genotype remains interpretable downstream.  VCF permits FORMAT
            # fields to vary by record.
            if len(raw_cols) == 1:
                vcf_format = "GT:AD:DP:LN"
            else:
                vcf_format = "GT:AD:DP:UC:UV:LN"
            sample_columns = [self._collapse_caller_blocks(event, raw_cols)]

        sample = "\t".join(sample_columns)

        return (
            f"{chrom}\t{pos}\t{record_id}\t{ref}\t{alt}\t{qual}\t"
            f"{filter_status}\t{info}\t{vcf_format}\t{sample}\n"
        )

    def _collapse_caller_blocks(self, event, raw_cols):
        """Collapse caller-mode evidence into one standard VCF sample column.

        A single evidence block is not a synthesis step and is preserved
        directly.  Once multiple evidence blocks are present, the output GT is
        produced by the same order-independent caller consensus used by sample
        merge.  Caller AD values are not composable, so synthesized calls carry
        ``AD=.,.`` and ``DP=.``.  UC/UV expose the unique carrier callers and
        valid caller votes that produced the synthesized call.  LN is taken
        from the merged event rather than from an arbitrarily selected caller
        block.
        """
        fmt = (
            event.format
            if getattr(event, "format", None)
            else "GT:AD:LN:ST:QV:TY:ID:SC:REF:ALT:CO"
        )

        if len(raw_cols) == 1:
            return self._convert_svcf_sample_to_vcf(raw_cols[0], fmt)

        consensus = resolve_caller_svcf_consensus(fmt, raw_cols, event.info)
        ln = self._event_length_for_synthesized_call(event)
        return (
            f"{consensus.gt}:.,.:.:"
            f"{consensus.unique_carrier_callers}:"
            f"{consensus.valid_caller_votes}:{ln}"
        )

    @staticmethod
    def _event_length_for_synthesized_call(event):
        """Return merged-event length for a synthesized caller-mode call."""
        svlen = getattr(event, "info", {}).get("SVLEN")
        if svlen not in (None, "", ".", True):
            try:
                return str(abs(int(str(svlen).strip())))
            except (TypeError, ValueError):
                pass

        sv_type = getattr(event, "sv_type", "")
        if sv_type not in {"BND", "TRA"}:
            try:
                start = int(getattr(event, "pos"))
                end = int(getattr(event, "end_pos"))
                if end != start:
                    return str(abs(end - start))
            except (TypeError, ValueError, AttributeError):
                pass

        return "."

    @staticmethod
    def _normalized_sample_value(value, default):
        if value in (None, ""):
            return default
        return str(value)

    def _convert_svcf_sample_to_vcf(
        self,
        svcf_sample,
        format_field,
        *,
        include_consensus_fields=False,
    ):
        """Convert one SVCF block using FORMAT keys rather than positions."""
        parsed = parse_svcf_sample_block(format_field, svcf_sample)

        if include_consensus_fields:
            # SVCF 1.1 sample-mode layout placeholders use GT=0/0, UC=0,
            # UV=0 to preserve fixed-width columns.  UV=0 explicitly means
            # there was no valid caller vote, so downstream VCF must expose
            # that state as missing rather than homozygous reference.  Legacy
            # sample layouts without UV are intentionally left unchanged.
            gt = downstream_sample_gt(
                parsed,
                unobserved_sample_gt=self.unobserved_sample_gt,
            )
        else:
            gt = self._normalized_sample_value(parsed.get("GT"), "./.")

        ad = self._normalized_sample_value(parsed.get("AD"), ".,.")
        ln = self._normalized_sample_value(parsed.get("LN"), ".")
        dp = self._calculate_dp(ad)

        if include_consensus_fields:
            uc = self._normalized_sample_value(parsed.get("UC"), ".")
            uv = self._normalized_sample_value(parsed.get("UV"), ".")
            return f"{gt}:{ad}:{dp}:{uc}:{uv}:{ln}"

        return f"{gt}:{ad}:{dp}:{ln}"

    def _get_alt(self, event):
        if event.alt and event.alt not in ("N", "."):
            return event.alt
        return f"<{event.sv_type}>"

    def _calculate_dp(self, ad):
        if ad in (None, "", "."):
            return "."

        parts = str(ad).split(",")
        if not parts or any(
            part == "." or not part.strip().isdigit()
            for part in parts
        ):
            return "."

        return sum(int(part) for part in parts)
