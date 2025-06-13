import logging

from .base import Converter


class NonBNDConverter(Converter):
    """A converter for non-BND events that enriches the event with additional info."""

    def _safe_int(self, value, default=0):
        """Safely convert value to int, return default if conversion fails."""
        try:
            if value == "." or value is None:
                return default
            return int(value)
        except (ValueError, TypeError):
            return default

    def convert(self, event):
        """Simplifies SVTYPE to one of the specified types and sets CHR2 and SVMETHOD for non-BND events."""
        # Check and simplify SVTYPE
        svtype = event.info.get("SVTYPE", "")
        if ":" in svtype:
            # Only take the part before the colon as SVTYPE
            svtype = svtype.split(":")[0]

        # Ensure the simplified SVTYPE is one of the allowed five types
        if svtype in ["INV", "INS", "DEL", "TRA", "DUP"]:
            event.info["SVTYPE"] = svtype
        else:
            # If SVTYPE is not among the specified types, print a warning or keep it unchanged
            logging.warning(f"Simplified SVTYPE '{svtype}' is not one of the allowed types. Keeping original SVTYPE.")

        # Only set CHR2 if it does not already exist
        if "CHR2" not in event.info:
            event.info["CHR2"] = event.chrom

        event.info["SVMETHOD"] = "OctopuSV"

        # Calculate SVLEN if not present and event is not TRA
        if "SVLEN" not in event.info and svtype in ["DEL", "DUP", "INV"]:
            try:
                pos = self._safe_int(event.pos)
                end = self._safe_int(event.info.get("END"))
                if pos != 0 and end != 0:  # Only calculate if we have valid positions
                    event.info["SVLEN"] = str(abs(end - pos))
                else:
                    event.info["SVLEN"] = "."
            except (ValueError, KeyError):
                logging.warning(f"Unable to calculate SVLEN for event {event.id}")
                event.info["SVLEN"] = "."

        # Set SVLEN to "." for TRA events
        if svtype == "TRA":
            event.info["SVLEN"] = "."

        # Deal with the "END" problem in INS
        if svtype == "INS":
            pos = self._safe_int(event.pos)
            svlen = self._safe_int(event.info.get("SVLEN", "0"))
            # If we have a valid position and SVLEN, calculate END
            if pos != 0 and svlen != 0:
                event.info["END"] = str(pos + svlen)
            else:
                # If we can't calculate a valid END position, use pos as END
                event.info["END"] = str(pos)

        # Ensure SVLEN exists in all cases
        if "SVLEN" not in event.info:
            event.info["SVLEN"] = "."
