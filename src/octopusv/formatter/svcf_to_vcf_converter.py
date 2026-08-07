import logging
import re

from octopusv.utils.genotype_resolver import resolve_multi_caller_genotype

logging.basicConfig(level=logging.INFO)


class SVCFtoVCFConverter:
    """Convert OctopuSV SVCF records back to VCF4.2-compatible records.

    Design goals:
      - Preserve standard VCF compatibility.
      - Preserve OctopuSV provenance in INFO, especially SOURCES and SOURCE_IDS.
      - Keep sample-mode sample columns intact.
      - Collapse caller-mode evidence blocks into one standard VCF sample column.
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

    def __init__(self, events, input_svcf_file):
        self.events = events
        self.input_svcf_file = input_svcf_file
        self.mode, self.sample_names = self._detect_mode_and_samples()

    @staticmethod
    def _header_id(line: str, kind: str) -> str | None:
        """Extract the ID from a VCF-style header definition line.

        Examples:
            ##INFO=<ID=SVTYPE,...> -> SVTYPE
            ##FORMAT=<ID=GT,...>   -> GT
            ##FILTER=<ID=PASS,...> -> PASS
        """
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
        """Make an INFO value safe for VCF text output.

        OctopuSV provenance values normally contain safe characters such as
        letters, digits, '.', '_', '-', ':', and commas. We keep commas because
        Number=. INFO values use comma-separated lists. Whitespace, semicolons,
        and equal signs would break INFO parsing, so we replace them.
        """
        text = str(value)
        text = re.sub(r"\s+", "_", text)
        text = text.replace(";", ",")
        text = text.replace("=", "_")
        return text

    def _detect_mode_and_samples(self):
        """Detect sample names from the #CHROM header.

        Multi vs single is determined by the number of sample columns in the
        #CHROM line. The returned `mode` is kept for backward compatibility,
        but sample output is driven by `self.sample_names`.
        """
        sample_names = ["Sample"]

        with open(self.input_svcf_file, encoding="utf-8-sig") as f:
            for line in f:
                if line.startswith("#CHROM"):
                    parts = line.strip().split("\t")
                    if len(parts) > 9:
                        sample_names = parts[9:]
                    break

        mode = "multi" if len(sample_names) > 1 else "single"
        return mode, sample_names

    def convert(self):
        vcf_content = self._generate_vcf_header()
        for event in self.events:
            vcf_content += self._convert_event_to_vcf(event)
        return vcf_content

    def _extract_original_definitions_from_svcf(self):
        """Extract non-OctopuSV-default header definitions from the SVCF file.

        We keep user/original caller definitions where possible, but skip
        OctopuSV default INFO/FORMAT/FILTER definitions because this converter
        writes a clean VCF4.2 header for the converted output.
        """
        header_definitions = {
            "filter_lines": [],
            "info_lines": [],
            "format_lines": [],
        }

        with open(self.input_svcf_file, encoding="utf-8-sig") as f:
            for line in f:
                line = line.strip()
                if line.startswith("#CHROM") or not line.startswith("##"):
                    break

                if line.startswith("##FILTER="):
                    filter_id = self._header_id(line, "FILTER")
                    if filter_id != "PASS":
                        header_definitions["filter_lines"].append(line)

                elif line.startswith("##INFO="):
                    info_id = self._header_id(line, "INFO")
                    if info_id is None or info_id not in self.DEFAULT_INFO_IDS:
                        header_definitions["info_lines"].append(line)

                elif line.startswith("##FORMAT="):
                    format_id = self._header_id(line, "FORMAT")
                    if format_id is None or format_id not in self.DEFAULT_FORMAT_IDS:
                        header_definitions["format_lines"].append(line)

        return header_definitions

    def _read_contig_lines(self) -> list[str]:
        """Read ##contig lines from the input SVCF header."""
        contig_lines = []
        with open(self.input_svcf_file, encoding="utf-8-sig") as f:
            for line in f:
                if line.startswith("##contig"):
                    contig_lines.append(line.rstrip("\n"))
                elif line.startswith("#CHROM"):
                    break
        return contig_lines

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

        # FILTER definitions.
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

        # INFO definitions.
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

        # ALT definitions.
        existing_alt_ids: set[str] = set()
        header = self._append_unique_header_lines(
            header,
            self.STANDARD_ALT_LINES,
            "ALT",
            existing_alt_ids,
        )

        # FORMAT definitions.
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

        # REF is required by VCF. Keep N as a valid reference placeholder rather
        # than converting it to '.', which is not valid for the REF column.
        ref = event.ref if event.ref and event.ref != "." else "N"

        alt = self._get_alt(event)
        qual = event.quality if hasattr(event, "quality") else "."
        filter_status = event.filter if hasattr(event, "filter") else "PASS"

        # Build standard + provenance-preserving INFO.
        info_fields = [f"SVTYPE={event.sv_type}"]

        # END handling.
        if event.sv_type in {"BND", "TRA"}:
            # BND/TRA mate coordinates are encoded in ALT.
            # Do not emit cross-chromosomal mate positions as INFO/END.
            pass
        elif event.sv_type == "INS":
            # VCF-compatible insertion representation.
            # SVCF may use END=POS+SVLEN internally.
            info_fields.append(f"END={event.pos}")
        else:
            if hasattr(event, "end_pos") and event.end_pos is not None:
                info_fields.append(f"END={event.end_pos}")

        # SVLEN handling.
        svlen_value = event.info.get("SVLEN")
        if not self._is_missing_info_value(svlen_value):
            try:
                svlen = int(str(svlen_value).strip())
                if event.sv_type == "DEL":
                    svlen = -abs(svlen)
                elif event.sv_type == "INS":
                    svlen = abs(svlen)

                # BND does not have a reliable linear SVLEN.
                if event.sv_type != "BND":
                    info_fields.append(f"SVLEN={svlen}")
            except ValueError:
                # Keep conversion VCF-compatible by skipping invalid SVLEN.
                pass

        # Preserve SVCF/OctopuSV INFO fields that remain meaningful in VCF.
        # SOURCES and SOURCE_IDS are essential provenance fields. They explain
        # which callers/samples supported a merged record, especially after
        # caller-mode evidence blocks are collapsed into one VCF sample column.
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

        # Resolve sample column(s).
        #
        # sample mode:
        #   Header has >1 sample. Emit one VCF column per sample, preserving all
        #   sample columns.
        #
        # caller mode:
        #   Header has one SAMPLE column, but each record may have multiple
        #   caller evidence blocks. Standard VCF cannot store multiple caller
        #   blocks as separate sample columns, so we collapse to one SAMPLE
        #   genotype and preserve caller provenance in INFO/SOURCES/SOURCE_IDS.
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
            sample_columns = [self._convert_svcf_sample_to_vcf(s) for s in raw_cols]
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
        """Collapse caller-mode evidence blocks into one VCF sample column.

        Uses the shared multi-caller voting rule to pick the representative
        genotype, then takes AD/LN from the first block carrying that genotype.
        This keeps svcf2vcf and stat consistent.
        """
        if len(raw_cols) == 1:
            return self._convert_svcf_sample_to_vcf(raw_cols[0])

        fmt = (
            event.format
            if getattr(event, "format", None)
            else "GT:AD:LN:ST:QV:TY:ID:SC:REF:ALT:CO"
        )
        info_str = ";".join(
            f"{k}={v}" if v is not True else k
            for k, v in event.info.items()
        )
        winning_gt = resolve_multi_caller_genotype(fmt, raw_cols, info_str)

        if winning_gt is not None:
            for col in raw_cols:
                parts = col.split(":")
                if parts and parts[0] == winning_gt:
                    return self._convert_svcf_sample_to_vcf(col)

        return self._convert_svcf_sample_to_vcf(raw_cols[0])

    def _convert_svcf_sample_to_vcf(self, svcf_sample):
        """Convert one SVCF sample/evidence block to GT:AD:DP:LN.

        SVCF FORMAT is usually:
            GT:AD:LN:ST:QV:TY:ID:SC:REF:ALT:CO

        Only GT, AD, and LN are needed here. They are the first three fields, so
        a simple split is safe even if later fields such as ID or ALT contain ':'.
        """
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