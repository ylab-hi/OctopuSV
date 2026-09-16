class SVCFtoBEDConverter:
    def __init__(self, events, minimal=False):
        self.events = events
        self.minimal = minimal

    def convert(self):
        bed_content = "" if self.minimal else 'track name=SVs description="Structural Variants from SVCF"\n'
        for event in self.events:
            bed_content += self._convert_event_to_bed(event)
        return bed_content

    def _convert_event_to_bed(self, event):
        """Convert one SVCF event to BED, failing loudly if it is not representable."""
        try:
            chrom = event.chrom
            start = int(event.pos) - 1

            if event.sv_type == "INS":
                end = start + 1
                svlen = event.info.get("SVLEN", "0")
            elif event.sv_type in ("TRA", "BND"):
                # BED cannot represent the mate breakpoint on another contig.
                # Represent the local breakpoint as a 1 bp interval instead of
                # inventing a same-contig span from CHR2/END.
                end = start + 1
                chr2 = event.info.get("CHR2", ".")
                end2 = event.info.get("END", ".")
                svlen = "."
                name = f"{event.sv_id}_{event.sv_type}_{chr2}:{end2}"
            else:
                end = int(event.info.get("END", start + 1))
                svlen = event.info.get("SVLEN", str(end - start))

            if self.minimal:
                return f"{chrom}\t{start}\t{end}\n"

            if event.sv_type not in ("TRA", "BND"):
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
        except (AttributeError, TypeError, ValueError) as exc:
            event_id = getattr(event, "sv_id", "<unknown>")
            raise ValueError(
                f"Could not convert SVCF event {event_id!r} to BED: {exc}"
            ) from exc
