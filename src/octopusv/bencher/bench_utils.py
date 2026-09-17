import json
from pathlib import Path

from octopusv.utils.text_io import open_text_auto
from octopusv.utils.vcf_info import format_vcf_info_item


_SAFE_HEADER_PREFIXES = (
    "##contig=",
    "##INFO=",
    "##FILTER=",
    "##ALT=",
    "##reference=",
)

_KNOWN_INFO_DEFINITIONS = {
    "SVTYPE": '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Structural variant type">',
    "END": '##INFO=<ID=END,Number=1,Type=Integer,Description="End position">',
    "SVLEN": '##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Structural variant length">',
    "CHR2": '##INFO=<ID=CHR2,Number=1,Type=String,Description="Second contig">',
    "SOURCES": '##INFO=<ID=SOURCES,Number=.,Type=String,Description="Supporting source labels">',
    "SOURCE_IDS": '##INFO=<ID=SOURCE_IDS,Number=.,Type=String,Description="Supporting source record IDs">',
}


def calculate_metrics(results: dict[str, list]):
    """Calculate benchmark metrics including precision, recall, and F1 score."""
    tp = len(results["tp_call"])
    fp = len(results["fp"])
    fn = len(results["fn"])

    precision = tp / (tp + fp) if tp + fp > 0 else 0
    recall = tp / (tp + fn) if tp + fn > 0 else 0
    f1 = 2 * (precision * recall) / (precision + recall) if precision + recall > 0 else 0

    return {"TP": tp, "FP": fp, "FN": fn, "Precision": precision, "Recall": recall, "F1": f1}


def read_safe_vcf_meta(path: Path | str) -> list[str]:
    """Return standard VCF metadata safe to inherit into benchmark outputs.

    SVCF identity and FORMAT declarations are deliberately omitted because
    benchmark result files contain only the eight fixed VCF columns and are
    not SVCF sample/evidence matrices.
    """
    fileformat = None
    inherited: list[str] = []
    seen: set[str] = set()

    with open_text_auto(path) as handle:
        for raw_line in handle:
            line = raw_line.rstrip("\r\n")
            if line.startswith("##fileformat="):
                fileformat = line
                continue
            if line.startswith(_SAFE_HEADER_PREFIXES):
                if line not in seen:
                    inherited.append(line)
                    seen.add(line)
                continue
            if line.startswith("#CHROM") or (line and not line.startswith("#")):
                break

    return [fileformat or "##fileformat=VCFv4.2", *inherited]


def _declared_info_ids(meta_lines: list[str]) -> set[str]:
    ids: set[str] = set()
    for line in meta_lines:
        if not line.startswith("##INFO=<ID="):
            continue
        rest = line[len("##INFO=<ID=") :]
        info_id = rest.split(",", 1)[0].split(">", 1)[0]
        if info_id:
            ids.add(info_id)
    return ids


def _event_info_items(event) -> dict:
    if isinstance(event, tuple):
        _, _, end_chrom, _, end_pos, source_file, _ = event
        return {
            "SVTYPE": "TRA",
            "END": end_pos,
            "CHR2": end_chrom,
            "SOURCES": source_file,
        }
    return dict(getattr(event, "info", {}) or {})


def _infer_info_definition(key: str, values: list[object]) -> str:
    if key in _KNOWN_INFO_DEFINITIONS:
        return _KNOWN_INFO_DEFINITIONS[key]

    nonmissing = [value for value in values if value not in (None, "", ".")]
    if nonmissing and all(value is True for value in nonmissing):
        return f'##INFO=<ID={key},Number=0,Type=Flag,Description="Benchmark-preserved INFO field">'

    number = "." if any("," in str(value) for value in nonmissing) else "1"
    info_type = "String"
    if nonmissing:
        try:
            for value in nonmissing:
                int(str(value))
            info_type = "Integer"
        except (TypeError, ValueError):
            try:
                for value in nonmissing:
                    float(str(value))
                info_type = "Float"
            except (TypeError, ValueError):
                info_type = "String"

    return (
        f'##INFO=<ID={key},Number={number},Type={info_type},'
        'Description="Benchmark-preserved INFO field">'
    )


def _complete_info_meta(meta_lines: list[str], events: list[tuple | object]) -> list[str]:
    """Ensure every INFO key emitted by benchmark has a header definition."""
    result = list(meta_lines)
    declared = _declared_info_ids(result)
    values_by_key: dict[str, list[object]] = {}
    for event in events:
        for key, value in _event_info_items(event).items():
            values_by_key.setdefault(str(key), []).append(value)

    for key in sorted(values_by_key):
        if key in declared:
            continue
        result.append(_infer_info_definition(key, values_by_key[key]))
        declared.add(key)

    return result


def write_vcf(
    file_path: Path,
    events: list[tuple | object],
    *,
    source_meta_lines: list[str] | None = None,
):
    """Write benchmark events as a self-describing, sample-free VCF."""
    meta_lines = _complete_info_meta(
        list(source_meta_lines or ["##fileformat=VCFv4.2"]),
        events,
    )

    with file_path.open("w") as f:
        for line in meta_lines:
            f.write(line.rstrip("\r\n") + "\n")
        f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        for event in events:
            if isinstance(event, tuple):  # TRA event
                chrom, start_chrom, end_chrom, start_pos, end_pos, source_file, bnd_pattern = event
                f.write(
                    f"{start_chrom}\t{start_pos}\t.\tN\t{bnd_pattern}\t.\tPASS\t"
                    f"SVTYPE=TRA;END={end_pos};CHR2={end_chrom};SOURCES={source_file}\n"
                )
            else:  # Other SV events
                f.write(
                    f"{event.chrom}\t{event.pos}\t{event.sv_id}\t{event.ref}\t{event.alt}\t"
                    f"{event.quality}\t{event.filter}\t{';'.join(format_vcf_info_item(k, v) for k, v in event.info.items())}\n"
                )


def write_summary(file_path: Path, metrics: dict):
    """Write benchmark summary metrics to JSON file."""
    with file_path.open("w") as f:
        json.dump(metrics, f, indent=2)
