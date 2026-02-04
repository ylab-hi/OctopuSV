class SVCFtoVCFConverter:
    def __init__(self, events, input_svcf_file):
        self.events = events
        self.input_svcf_file = input_svcf_file

    def convert(self):
        vcf_content = self._generate_vcf_header()
        for event in self.events:
            vcf_content += self._convert_event_to_vcf(event)
        return vcf_content

    def _extract_original_definitions_from_svcf(self):
        """Extract original header definitions from SVCF file"""
        header_definitions = {"filter_lines": [], "info_lines": [], "format_lines": [], "other_lines": []}

        with open(self.input_svcf_file) as f:
            for line in f:
                line = line.strip()
                if line.startswith("#CHROM") or not line.startswith("##"):
                    break

                if line.startswith("##FILTER=") and "PASS" not in line:
                    # Skip OctopuSV default PASS filter, keep original filters
                    header_definitions["filter_lines"].append(line)
                elif line.startswith("##INFO=") and not any(
                    x in line
                    for x in ["SVTYPE", "END", "SVLEN", "CHR2", "SUPPORT", "SVMETHOD", "STRAND", "RNAMES", "RTID", "AF"]
                ):
                    # Keep original INFO fields, skip OctopuSV defaults
                    header_definitions["info_lines"].append(line)
                elif line.startswith("##FORMAT=") and not any(
                    x in line for x in ["GT", "AD", "DP", "LN", "ST", "QV", "TY", "ID", "SC", "REF", "ALT", "CO"]
                ):
                    # Keep original FORMAT fields, skip our defaults
                    header_definitions["format_lines"].append(line)

        return header_definitions

    def _generate_vcf_header(self):
        """Generate VCF header with original definitions from SVCF"""
        # Extract original definitions from SVCF
        original_defs = self._extract_original_definitions_from_svcf()

        # Read contig info from SVCF
        contig_lines = ""
        with open(self.input_svcf_file) as f:
            for line in f:
                if line.startswith("##contig"):
                    contig_lines += line
                elif line.startswith("#CHROM"):
                    break

        # Build header with original definitions
        header = f"""##fileformat=VCFv4.2
{contig_lines}"""

        # Add original FILTER definitions first
        for filter_line in original_defs["filter_lines"]:
            header += filter_line + "\n"

        # Add standard FILTER if not already present
        if not any("PASS" in line for line in original_defs["filter_lines"]):
            header += '##FILTER=<ID=PASS,Description="All filters passed">\n'

        # Add original INFO definitions
        for info_line in original_defs["info_lines"]:
            header += info_line + "\n"

        # Add our standard INFO definitions
        header += """##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">
##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant described in this record">
##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Difference in length between REF and ALT alleles">
##INFO=<ID=CHR2,Number=1,Type=String,Description="Chromosome for end coordinate">
##INFO=<ID=SUPPORT,Number=1,Type=Integer,Description="Number of reads supporting this variant">
##INFO=<ID=SVMETHOD,Number=1,Type=String,Description="Method used to detect SV">
##INFO=<ID=STRAND,Number=1,Type=String,Description="Strand orientation of the SV">
##ALT=<ID=DEL,Description="Deletion">
##ALT=<ID=INV,Description="Inversion">
##ALT=<ID=DUP,Description="Duplication">
##ALT=<ID=INS,Description="Insertion">
##ALT=<ID=TRA,Description="Translocation">
##ALT=<ID=BND,Description="Breakend">
"""

        # Add original FORMAT definitions
        for format_line in original_defs["format_lines"]:
            header += format_line + "\n"

        # Add our standard FORMAT definitions
        header += """##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth">
##FORMAT=<ID=LN,Number=1,Type=Integer,Description="Length of SV">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSample\n"""

        return header

    def _convert_event_to_vcf(self, event):
        chrom = event.chrom
        pos = event.pos
        id = event.sv_id
        ref = event.ref if event.ref != "N" else "."
        alt = self._get_alt(event)
        qual = event.quality if hasattr(event, "quality") else "."
        filter = event.filter if hasattr(event, "filter") else "PASS"

        # Start with SVTYPE
        info_fields = [f"SVTYPE={event.sv_type}"]

        # Handle END field - avoid END < POS coordinate issues
        if event.sv_type in ["BND"]:
            # BND events don't use END field, coordinates are in ALT field
            pass
        elif event.sv_type == "INS":
            # For insertions, END should equal POS
            end_pos = event.pos
            info_fields.append(f"END={end_pos}")
        elif event.sv_type in ["TRA"]:
            # For translocation events, always include END if it exists
            if hasattr(event, "end_pos") and event.end_pos is not None:
                info_fields.append(f"END={event.end_pos}")
        else:
            # For DEL, DUP, INV - use actual end position
            end_pos = event.end_pos
            info_fields.append(f"END={end_pos}")

        # Handle SVLEN calculation
        if "SVLEN" in event.info and event.info["SVLEN"] != ".":
            svlen_str = event.info["SVLEN"]
            try:
                svlen = int(svlen_str.strip())
                if event.sv_type == "DEL":
                    svlen = -abs(svlen)
                elif event.sv_type == "INS":
                    svlen = abs(svlen)
                elif event.sv_type == "BND":
                    # BND events typically don't have meaningful SVLEN
                    pass  # Don't add SVLEN for BND
                else:
                    # DUP, INV, TRA keep original value
                    pass

                # Only add SVLEN for non-BND events
                if event.sv_type != "BND":
                    info_fields.append(f"SVLEN={svlen}")
            except ValueError:
                pass  # Skip invalid SVLEN values

        # Add other INFO fields only if they have meaningful values
        if "CHR2" in event.info and event.info["CHR2"] != "." and event.info["CHR2"]:
            info_fields.append(f"CHR2={event.info['CHR2']}")
        if "SUPPORT" in event.info and event.info["SUPPORT"] != "." and event.info["SUPPORT"]:
            info_fields.append(f"SUPPORT={event.info['SUPPORT']}")
        if "SVMETHOD" in event.info and event.info["SVMETHOD"] != "." and event.info["SVMETHOD"]:
            info_fields.append(f"SVMETHOD={event.info['SVMETHOD']}")
        if "STRAND" in event.info and event.info["STRAND"] != "." and event.info["STRAND"]:
            info_fields.append(f"STRAND={event.info['STRAND']}")

        info = ";".join(info_fields)

        format = "GT:AD:DP:LN"
        gt = event.sample.get("GT", "./.")
        ad = event.sample.get("AD", ".,.")
        dp = self._calculate_dp(ad)
        ln = event.sample.get("LN", ".")

        sample = f"{gt}:{ad}:{dp}:{ln}"

        return f"{chrom}\t{pos}\t{id}\t{ref}\t{alt}\t{qual}\t{filter}\t{info}\t{format}\t{sample}\n"

    def _get_alt(self, event):
        if event.alt and event.alt not in ("N", "."):
            return event.alt
        return f"<{event.sv_type}>"

    def _calculate_dp(self, ad):
        try:
            return sum(int(x) for x in ad.split(",") if x != "." and x.strip().isdigit())
        except ValueError:
            return "."
