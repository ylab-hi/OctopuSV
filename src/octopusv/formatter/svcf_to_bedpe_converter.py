class SVCFtoBEDPEConverter:
    def __init__(self, events, minimal=False):
        self.events = events
        self.minimal = minimal

    def convert(self):
        if self.minimal:
            header = "#chrom1\tstart1\tend1\tchrom2\tstart2\tend2\n"
        else:
            header = "#chrom1\tstart1\tend1\tchrom2\tstart2\tend2\tname\tscore\tstrand1\tstrand2\tsvtype\tsvlen\n"

        bedpe_content = header
        for event in self.events:
            bedpe_content += self._convert_event_to_bedpe(event)
        return bedpe_content

    def _convert_event_to_bedpe(self, event):
        """Convert one SVCF event to BEDPE, failing loudly if not representable."""
        try:
            chrom1 = event.chrom
            start1 = int(event.pos) - 1
            end1 = int(event.pos)

            if event.sv_type in ("TRA", "BND"):
                chrom2 = event.info.get("CHR2")
                end_value = event.info.get("END")
                if chrom2 in (None, "", ".") or end_value in (None, "", "."):
                    raise ValueError(
                        f"{event.sv_type} requires explicit CHR2 and numeric END "
                        "for BEDPE conversion"
                    )
            else:
                chrom2 = event.info.get("CHR2", chrom1)
                end_value = event.info.get("END", end1)

            start2 = int(end_value) - 1
            end2 = int(end_value)

            if self.minimal:
                return f"{chrom1}\t{start1}\t{end1}\t{chrom2}\t{start2}\t{end2}\n"

            name = f"{event.sv_id}_{event.sv_type}"
            score = event.info.get("SUPPORT", event.quality if hasattr(event, "quality") else "1")

            strand1 = "+"
            strand2 = "-" if event.sv_type in ["INV", "TRA"] else "+"
            strand_info = event.info.get("STRAND", "")
            if strand_info == "+-":
                strand1, strand2 = "+", "-"
            elif strand_info == "-+":
                strand1, strand2 = "-", "+"
            elif strand_info == "++":
                strand1, strand2 = "+", "+"
            elif strand_info == "--":
                strand1, strand2 = "-", "-"

            svtype = event.sv_type
            svlen = event.info.get("SVLEN", ".")
            return f"{chrom1}\t{start1}\t{end1}\t{chrom2}\t{start2}\t{end2}\t{name}\t{score}\t{strand1}\t{strand2}\t{svtype}\t{svlen}\n"
        except (AttributeError, TypeError, ValueError) as exc:
            event_id = getattr(event, "sv_id", "<unknown>")
            raise ValueError(
                f"Could not convert SVCF event {event_id!r} to BEDPE: {exc}"
            ) from exc
