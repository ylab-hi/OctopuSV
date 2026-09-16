"""Single-pass SVCF reader shared by all stat analyzers.

The old stat code opened the file once per analyzer (5x IO on large files) and
each analyzer parsed INFO its own (fragile) way. This module reads the file
ONCE, parsing every record into a small structured dict, and every analyzer
consumes that same in-memory list. INFO parsing is centralized and safe
(handles flags and values containing '=').
"""


def parse_info(info_str):
    """Parse an INFO column into a dict. Flags (no '=') map to True.

    Safe against values that themselves contain '=' (splits on the first only),
    unlike dict(item.split('=')).
    """
    info = {}
    for item in info_str.split(";"):
        if not item:
            continue
        if "=" in item:
            key, value = item.split("=", 1)
            info[key] = value
        else:
            info[item] = True
    return info


class SVRecord:
    """A minimal structured view of one SVCF data row.

    Only the fields the analyzers need are kept. Sample columns are retained
    raw (one string per caller/sample) so genotype resolution can run later.
    """

    __slots__ = ("chrom", "pos", "sv_id", "qual", "filter", "info",
                 "format", "sample_cols", "svtype")

    def __init__(self, fields):
        self.chrom = fields[0]
        self.pos = fields[1]
        self.sv_id = fields[2]
        self.qual = fields[5]
        self.filter = fields[6]
        self.info = parse_info(fields[7])
        self.format = fields[8] if len(fields) > 8 else ""
        self.sample_cols = fields[9:] if len(fields) > 9 else []
        self.svtype = self.info.get("SVTYPE", "Unknown")


def read_records(input_file):
    """Read an SVCF file once. Returns (records, sample_names, mode).

    ``mode`` is ``sample`` when ``##OctopuSV_mode=multi`` is declared.
    Header-only legacy multi-sample files without the marker are still treated
    as sample mode when they contain more than one trailing column. Otherwise
    the file is caller/single mode.
    """
    records = []
    sample_names = []
    has_multi_marker = False
    with open(input_file) as fh:
        for line in fh:
            if line.rstrip("\r\n") == "##OctopuSV_mode=multi":
                has_multi_marker = True
                continue
            if line.startswith("#CHROM"):
                header = line.rstrip("\n").split("\t")
                sample_names = header[9:] if len(header) > 9 else []
                continue
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 8:
                continue
            records.append(SVRecord(fields))

    mode = "sample" if has_multi_marker or len(sample_names) > 1 else "caller"
    return records, sample_names, mode
