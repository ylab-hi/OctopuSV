from enum import Enum

from octopusv.utils.construct_sample_string import construct_sample_string


class SVType(Enum):
    """A class to represent the SV type."""

    DEL = "DEL"
    DUP = "DUP"
    INV = "INV"
    INS = "INS"
    TRA = "TRA"
    BND = "BND"


class SVEvent:
    """represent each SV event."""

    def __init__(self, chrom, pos, id, ref, alt, qual, filter, info, format="GT", sample="0/1"):
        self.chrom = chrom
        self.pos = int(pos)
        self.id = id
        self.ref = ref
        self.alt = alt
        self.orig_alt = alt
        self.qual = qual
        self.filter = filter
        self.info = self._parse_info(info)  # self.info will be a dict
        self.format = format
        self.sample = sample

    def _parse_info(self, info):  # Parse the info field according to your requirement
        info_dict = {}
        for item in info.split(";"):
            parts = item.split("=")
            if len(parts) == 2:  # Just in case that the INFO doesn't contain = symbol
                key, value = parts
                if key.lower() == "strand":
                    key = "STRAND"
                info_dict[key] = value
            else:
                info_dict[parts[0]] = None

        for support_key in ["SUPP", "SUPPREAD", "WEIGHT", "RE"]:
            if support_key in info_dict:
                info_dict["SUPPORT"] = info_dict[support_key]
                break

        # Handle different field names for RNAMES (e.g., SVIM uses READS)
        for rnames_key in ["READS", "RNAMES"]:
            if rnames_key in info_dict:
                info_dict["RNAMES"] = info_dict[rnames_key]
                break

        # If RNAMES not in info, set it to "."
        if "RNAMES" not in info_dict:
            info_dict["RNAMES"] = "."

        return info_dict

    def __getitem__(self, key):
        # enable you to access object as dict
        return self.info[key]

    def __setitem__(self, key, value):
        self.info[key] = value

    def is_duplication(self):
        # Whether is duplication
        return self.info["SVTYPE"] == "DUP"

    def is_inversion(self):
        # Whether is inversion
        return self.info["SVTYPE"] == "INV"

    def is_insertion(self):
        #  Whether is insertion
        return self.info["SVTYPE"] == "INS"

    def is_TRA(self):
        #  Whether is TRA
        return self.info["SVTYPE"] == "TRA"

    def is_BND(self):
        #  Whether is BND
        return self.info["SVTYPE"] == "BND"

    def __str__(self):
        """Convert the SVEvent to a string in SVCF format."""
        # Process 'PR' for 'SUPPORT' if not already set and if FORMAT/SAMPLE fields exist
        if "SUPPORT" not in self.info and self.format != "GT":
            format_parts = self.format.split(":")
            sample_parts = self.sample.split(":")
            format_sample_dict = dict(zip(format_parts, sample_parts, strict=False))
            # Handling both 'PR' only and 'PR:SR' cases.
            pr_key = "PR" if "PR" in format_sample_dict else None
            if pr_key:
                pr_values = format_sample_dict[pr_key].split(",")
                if len(pr_values) > 1:
                    self.info["SUPPORT"] = pr_values[1]  # Use the alt allele's paired-read count

        # Fixed order for INFO fields, using '.' as a placeholder for missing values
        info_order = ["SVTYPE", "END", "SVLEN", "CHR2", "SUPPORT", "SVMETHOD", "RTID", "AF", "STRAND", "RNAMES"]
        info_str_parts = []
        for key in info_order:
            value = self.info.get(key, ".")
            if key == "SVLEN" and value != ".":  # Check if SVLEN is present and not placeholder
                try:
                    # Ensure the value is an integer and take its absolute value
                    value = str(abs(int(value)))
                except ValueError:
                    # In case value is not an integer, keep the original value
                    pass
            info_str_parts.append(f"{key}={value}")
        info_str = ";".join(info_str_parts)

        # Always include FORMAT and SAMPLE fields
        format_order = ["GT", "AD", "LN", "ST", "QV", "TY", "ID", "SC", "REF", "ALT", "CO"]
        format_str = ":".join(format_order)
        sample_str = construct_sample_string(self)

        return f"{self.chrom}\t{self.pos}\t{self.id}\t{self.ref}\t{self.alt}\t{self.qual}\t{self.filter}\t{info_str}\t{format_str}\t{sample_str}"
