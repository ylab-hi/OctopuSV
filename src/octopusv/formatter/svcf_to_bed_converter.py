class SVCFtoBEDConverter:
    def __init__(self, events, minimal=False):
        self.events = events
        self.minimal = minimal

    def convert(self):
        # Minimal format doesn't need track line
        bed_content = "" if self.minimal else 'track name=SVs description="Structural Variants from SVCF"\n'

        for event in self.events:
            bed_line = self._convert_event_to_bed(event)
            if bed_line:  # Only add valid conversions
                bed_content += bed_line
        return bed_content

    def _convert_event_to_bed(self, event):
        """Convert SVCF event to BED format Minimal BED format: chrom start end Standard BED format: chrom start end name score strand."""
        try:
            chrom = event.chrom
            start = int(event.pos) - 1  # BED is 0-based

            # Handle different SV types
            if event.sv_type == "INS":
                # For insertions, use a 1bp interval at insertion point
                end = start + 1
                svlen = event.info.get("SVLEN", "0")
            elif event.sv_type == "TRA":
                # For translocations, show both source and destination
                end = start + 1
                chr2 = event.info.get("CHR2", "")
                end2 = event.info.get("END", "")
                svlen = "."
                name = f"{event.sv_id}_TRA_{chr2}:{end2}"
            else:
                # For DEL, DUP, INV use the full interval
                end = int(event.info.get("END", start + 1))
                svlen = event.info.get("SVLEN", str(end - start))

            # Minimal format only includes first three columns
            if self.minimal:
                return f"{chrom}\t{start}\t{end}\n"

            # Full format includes all columns
            if event.sv_type != "TRA":
                name = f"{event.sv_id}_{event.sv_type}_{svlen}bp"

            score = event.info.get("SUPPORT", event.quality if hasattr(event, "quality") else "0")

            strand = event.info.get("STRAND", ".")
            if strand in ("+-", "-+"):
                strand = "-"
            elif strand == "++":
                strand = "+"
            elif strand not in ["+", "-"]:
                strand = "."

            return f"{chrom}\t{start}\t{end}\t{name}\t{score}\t{strand}\n"

        except (AttributeError, ValueError) as e:
            typer.echo(f"Warning: Could not convert event {event.sv_id}: {e!s}")
            return None
