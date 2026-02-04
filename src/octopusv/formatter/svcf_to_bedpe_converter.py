import logging


class SVCFtoBEDPEConverter:
    def __init__(self, events, minimal=False):
        self.events = events
        self.minimal = minimal

    def convert(self):
        if self.minimal:
            # Minimal format: only essential columns
            header = "#chrom1\tstart1\tend1\tchrom2\tstart2\tend2\n"
        else:
            # Full format: all columns
            header = "#chrom1\tstart1\tend1\tchrom2\tstart2\tend2\tname\tscore\tstrand1\tstrand2\tsvtype\tsvlen\n"

        bedpe_content = header
        for event in self.events:
            bedpe_line = self._convert_event_to_bedpe(event)
            if bedpe_line:  # Only add valid conversions
                bedpe_content += bedpe_line
        return bedpe_content

    def _convert_event_to_bedpe(self, event):
        """Convert SVCF event to BEDPE format."""
        try:
            # Process first breakpoint
            chrom1 = event.chrom
            start1 = int(event.pos) - 1  # BEDPE is 0-based
            end1 = int(event.pos)

            # Process second breakpoint
            chrom2 = event.info.get("CHR2", chrom1)
            start2 = int(event.info.get("END", end1)) - 1
            end2 = int(event.info.get("END", end1))

            if self.minimal:
                # Minimal format: only coordinates
                return f"{chrom1}\t{start1}\t{end1}\t{chrom2}\t{start2}\t{end2}\n"

            # Full format: include all information
            name = f"{event.sv_id}_{event.sv_type}"
            score = event.info.get("SUPPORT", event.quality if hasattr(event, "quality") else "1")

            # Get strand information
            strand1 = "+"
            strand2 = "-" if event.sv_type in ["INV", "TRA"] else "+"

            strand_info = event.info.get("STRAND", "")
            if strand_info:
                if strand_info == "+-":
                    strand1, strand2 = "+", "-"
                elif strand_info == "-+":
                    strand1, strand2 = "-", "+"
                elif strand_info == "++":
                    strand1, strand2 = "+", "+"
                elif strand_info == "--":
                    strand1, strand2 = "-", "-"

            # Get SV type and length
            svtype = event.sv_type
            svlen = event.info.get("SVLEN", ".")

            return f"{chrom1}\t{start1}\t{end1}\t{chrom2}\t{start2}\t{end2}\t{name}\t{score}\t{strand1}\t{strand2}\t{svtype}\t{svlen}\n"

        except (AttributeError, ValueError) as e:
            logging.error(f"Warning: Could not convert event {event.sv_id}: {e!s}")
            return None
