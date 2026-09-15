import io
import logging
import os
import re
import uuid
from pathlib import Path

from octopusv.utils.genotype_resolver import (
    resolve_multi_caller_genotype,
    unique_source_segments,
)
from octopusv.utils.svcf_parser import SVCFEvent

logging.basicConfig(level=logging.INFO)


class SVCFtoVCFConverter:
    """Convert OctopuSV SVCF records back to VCF4.2-compatible records.

    Design goals:
      - Preserve standard VCF compatibility.
      - Preserve OctopuSV source/evidence metadata in INFO, especially
        SOURCES and SOURCE_IDS.
      - Keep sample-mode sample columns intact and in header order.
      - Collapse caller-mode evidence blocks into one standard VCF sample
        column using the shared unique-source genotype rule.
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

    def __init__(self, events=None, input_svcf_file=None):
        if input_svcf_file is None:
            raise ValueError("input_svcf_file is required")

        self.events = events
        self.input_svcf_file = str(input_svcf_file)

        (
            self._meta_lines,
            self._contig_lines,
            self._original_definitions,
            self.sample_names,
        ) = self._read_header_metadata()

        self.mode = "multi" if len(self.sample_names) > 1 else "single"

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

        Returns meta lines, contig lines, preserved original definitions, and
        the trailing #CHROM sample/evidence-column names.  This replaces the
        previous repeated header scans without changing header semantics.
        """
        meta_lines = []
        contig_lines = []
        original_definitions = {
            "filter_lines": [],
            "info_lines": [],
            "format_lines": [],
        }
        sample_names = ["Sample"]

        with open(self.input_svcf_file, encoding="utf-8-sig") as handle:
            for raw_line in handle:
                line = raw_line.rstrip("\r\n")

                if line.startswith("#CHROM"):
                    parts = line.split("\t")
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

                elif line.startswith("##FORMAT="):
                    format_id = self._header_id(line, "FORMAT")
                    if format_id is None or format_id not in self.DEFAULT_FORMAT_IDS:
                        original_definitions["format_lines"].append(line)

        return meta_lines, contig_lines, original_definitions, sample_names

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
        """Stream conversion to a temporary file, then atomically replace output.

        Streaming avoids holding the complete SVCF and VCF in memory.  The
        temporary-file/replace pattern preserves the previous all-or-nothing
        behavior: an error midway through conversion must not leave a partial
        VCF at the requested output path.

        The temporary file is created with normal file-creation permissions
        (subject to the process umask).  When replacing an existing output, its
        permission bits are preserved.
        """
        output_path = Path(output_file)
        output_dir = output_path.parent
        existing_mode = None
        if output_path.exists():
            existing_mode = output_path.stat().st_mode & 0o777

        tmp_path = None
        fd = None
        try:
            # Use O_EXCL so the temporary path cannot collide with an existing
            # file.  Mode 0o666 is masked by the caller's normal process umask,
            # matching ordinary text-file creation better than NamedTemporaryFile
            # (which forces 0o600).
            for _ in range(10):
                candidate = output_dir / (
                    f".{output_path.name}.{uuid.uuid4().hex}.tmp"
                )
                try:
                    fd = os.open(
                        candidate,
                        os.O_WRONLY | os.O_CREAT | os.O_EXCL,
                        0o666,
                    )
                    tmp_path = candidate
                    break
                except FileExistsError:
                    continue

            if fd is None or tmp_path is None:
                raise OSError(
                    f"Could not create temporary output beside {output_path}"
                )

            if existing_mode is not None:
                os.fchmod(fd, existing_mode)

            with os.fdopen(fd, mode="w", encoding="utf-8") as handle:
                fd = None
                self.write(handle)

            os.replace(tmp_path, output_path)

        except Exception:
            if fd is not None:
                os.close(fd)
            if tmp_path is not None:
                try:
                    tmp_path.unlink()
                except FileNotFoundError:
                    pass
            raise

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

        with open(self.input_svcf_file, encoding="utf-8-sig") as handle:
            for line in handle:
                if line.startswith("#"):
                    continue

                parts = line.rstrip("\r\n").split("\t")
                if len(parts) < 10:
                    # Preserve SVCFFileEventCreator's historical behavior for
                    # incomplete data rows: ignore them rather than inventing
                    # missing fixed/evidence columns.
                    continue

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

        existing_alt_ids: set[str] = set()
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
        chrom = event.chrom
        pos = event.pos
        record_id = event.sv_id

        ref = event.ref if event.ref and event.ref != "." else "N"

        alt = self._get_alt(event)
        qual = event.quality if hasattr(event, "quality") else "."
        filter_status = event.filter if hasattr(event, "filter") else "PASS"

        info_fields = [f"SVTYPE={event.sv_type}"]

        if event.sv_type in {"BND", "TRA"}:
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

        info = ";".join(info_fields) if info_fields else "."
        vcf_format = "GT:AD:DP:LN"

        raw_cols = getattr(event, "raw_sample_columns", None)
        if not raw_cols:
            gt = event.sample.get("GT", "./.")
            ad = event.sample.get("AD", ".,.")
            ln = event.sample.get("LN", ".")
            raw_cols = [f"{gt}:{ad}:{ln}"]

        expected = len(self.sample_names)

        if expected == 1:
            sample_columns = [self._collapse_caller_blocks(event, raw_cols)]
        else:
            sample_columns = [
                self._convert_svcf_sample_to_vcf(sample)
                for sample in raw_cols
            ]
            if len(sample_columns) != expected:
                raise ValueError(
                    f"Sample column count mismatch for {record_id} at {chrom}:{pos}: "
                    f"got {len(sample_columns)}, expected {expected} "
                    "(from #CHROM header)."
                )

        sample = "\t".join(sample_columns)

        return (
            f"{chrom}\t{pos}\t{record_id}\t{ref}\t{alt}\t{qual}\t"
            f"{filter_status}\t{info}\t{vcf_format}\t{sample}\n"
        )

    def _collapse_caller_blocks(self, event, raw_cols):
        """Collapse caller evidence blocks into one VCF sample column."""
        if len(raw_cols) == 1:
            return self._convert_svcf_sample_to_vcf(raw_cols[0])

        fmt = (
            event.format
            if getattr(event, "format", None)
            else "GT:AD:LN:ST:QV:TY:ID:SC:REF:ALT:CO"
        )
        info_str = ";".join(
            f"{key}={value}" if value is not True else key
            for key, value in event.info.items()
        )
        winning_gt = resolve_multi_caller_genotype(fmt, raw_cols, info_str)
        selected_blocks = unique_source_segments(raw_cols, info_str)

        # AD/LN must come from the same unique-source evidence set that
        # participated in genotype resolution.  Otherwise an excluded duplicate
        # block from the same caller can accidentally supply AD/LN merely because
        # it carries the winning GT.
        if winning_gt is not None:
            for _index, col in selected_blocks:
                parts = col.split(":")
                if parts and parts[0] == winning_gt:
                    return self._convert_svcf_sample_to_vcf(col)

        if selected_blocks:
            return self._convert_svcf_sample_to_vcf(selected_blocks[0][1])
        return self._convert_svcf_sample_to_vcf(raw_cols[0])

    def _convert_svcf_sample_to_vcf(self, svcf_sample):
        """Convert one SVCF sample/evidence block to GT:AD:DP:LN."""
        parts = svcf_sample.split(":")

        if len(parts) < 3:
            return "./.:.,.:0:."

        gt = parts[0] if parts[0] else "./."
        ad = parts[1] if len(parts) > 1 else ".,."
        ln = parts[2] if len(parts) > 2 else "."

        dp = self._calculate_dp(ad)

        return f"{gt}:{ad}:{dp}:{ln}"

    def _get_alt(self, event):
        if event.alt and event.alt not in ("N", "."):
            return event.alt
        return f"<{event.sv_type}>"

    def _calculate_dp(self, ad):
        try:
            return sum(
                int(x)
                for x in ad.split(",")
                if x != "." and x.strip().isdigit()
            )
        except ValueError:
            return "."
