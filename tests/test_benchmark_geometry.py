from __future__ import annotations

from pathlib import Path

from octopusv.bencher.sv_bencher import SVBencher
from octopusv.utils.svcf_schema import CALLER_FORMAT


def _sample_block(*, svtype: str, record_id: str, chrom: str, pos: int, end_chrom: str, end: int) -> str:
    length = abs(end - pos) if svtype not in {"TRA", "BND"} else "."
    return f"0/1:5,5:{length}:.:60:{svtype}:{record_id}:test:N:<{svtype}>:{chrom}_{pos}-{end_chrom}_{end}"


def _write_svcf(path: Path, records: list[dict]) -> None:
    contigs = sorted({record["chrom"] for record in records} | {record["chr2"] for record in records})
    lines = [
        "##fileformat=VCFv4.2",
        *[f"##contig=<ID={contig},length=5000000>" for contig in contigs],
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE",
    ]

    for record in records:
        svtype = record["svtype"]
        chrom = record["chrom"]
        pos = int(record["pos"])
        chr2 = record.get("chr2", chrom)
        end = int(record["end"])
        record_id = record["id"]
        filt = record.get("filter", "PASS")

        if svtype in {"TRA", "BND"}:
            alt = record.get("alt", f"N]{chr2}:{end}]")
            svlen = "."
        else:
            alt = record.get("alt", f"<{svtype}>")
            svlen = str(abs(end - pos))

        info = ";".join(
            [
                f"SVTYPE={svtype}",
                f"END={end}",
                f"SVLEN={svlen}",
                f"CHR2={chr2}",
                "SUPPORT=10",
                "SVMETHOD=test",
                "RTID=.",
                "AF=.",
                "STRAND=.",
                "RNAMES=.",
            ]
        )
        sample = _sample_block(
            svtype=svtype,
            record_id=record_id,
            chrom=chrom,
            pos=pos,
            end_chrom=chr2,
            end=end,
        )
        lines.append(
            f"{chrom}\t{pos}\t{record_id}\tN\t{alt}\t60\t{filt}\t{info}\t{CALLER_FORMAT}\t{sample}"
        )

    path.write_text("\n".join(lines) + "\n")


def _run(tmp_path: Path, truth: list[dict], calls: list[dict], **kwargs):
    truth_path = tmp_path / "truth.svcf"
    call_path = tmp_path / "call.svcf"
    output_dir = tmp_path / "bench"
    _write_svcf(truth_path, truth)
    _write_svcf(call_path, calls)

    bencher = SVBencher(truth_path, call_path, output_dir, **kwargs)
    bencher.run_benchmark()
    return bencher.results


def _counts(results):
    return (
        len(results["tp_call"]),
        len(results["fp"]),
        len(results["fn"]),
    )


def test_existing_interval_del_match_stays_tp(tmp_path):
    truth = [{"id": "t", "svtype": "DEL", "chrom": "chr1", "pos": 1000, "chr2": "chr1", "end": 1500}]
    calls = [{"id": "c", "svtype": "DEL", "chrom": "chr1", "pos": 1010, "chr2": "chr1", "end": 1510}]
    assert _counts(_run(tmp_path, truth, calls)) == (1, 0, 0)


def test_existing_interval_ins_match_stays_tp(tmp_path):
    truth = [{"id": "t", "svtype": "INS", "chrom": "chr1", "pos": 2000, "chr2": "chr1", "end": 2100}]
    calls = [{"id": "c", "svtype": "INS", "chrom": "chr1", "pos": 2010, "chr2": "chr1", "end": 2110}]
    assert _counts(_run(tmp_path, truth, calls)) == (1, 0, 0)


def test_interval_events_on_different_chromosomes_do_not_match(tmp_path):
    truth = [{"id": "t", "svtype": "DEL", "chrom": "chr1", "pos": 1000, "chr2": "chr1", "end": 1500}]
    calls = [{"id": "c", "svtype": "DEL", "chrom": "chr2", "pos": 1000, "chr2": "chr2", "end": 1500}]
    assert _counts(_run(tmp_path, truth, calls)) == (0, 1, 1)


def test_bnd_same_two_breakpoint_geometry_matches(tmp_path):
    truth = [{"id": "t", "svtype": "BND", "chrom": "chr1", "pos": 1000, "chr2": "chr2", "end": 2000}]
    calls = [{"id": "c", "svtype": "BND", "chrom": "chr1", "pos": 1010, "chr2": "chr2", "end": 2010}]
    assert _counts(_run(tmp_path, truth, calls)) == (1, 0, 0)


def test_bnd_wrong_mate_chromosome_does_not_match(tmp_path):
    truth = [{"id": "t", "svtype": "BND", "chrom": "chr1", "pos": 1000, "chr2": "chr2", "end": 2000}]
    calls = [{"id": "c", "svtype": "BND", "chrom": "chr1", "pos": 1000, "chr2": "chr3", "end": 2000}]
    assert _counts(_run(tmp_path, truth, calls)) == (0, 1, 1)


def test_bnd_second_breakpoint_outside_reference_distance_does_not_match(tmp_path):
    truth = [{"id": "t", "svtype": "BND", "chrom": "chr1", "pos": 1000, "chr2": "chr2", "end": 2000}]
    calls = [{"id": "c", "svtype": "BND", "chrom": "chr1", "pos": 1000, "chr2": "chr2", "end": 2601}]
    assert _counts(_run(tmp_path, truth, calls, reference_distance=500)) == (0, 1, 1)


def test_bnd_is_not_filtered_by_numeric_cross_chromosome_span(tmp_path):
    truth = [{"id": "t", "svtype": "BND", "chrom": "chr1", "pos": 100, "chr2": "chr2", "end": 1000000}]
    calls = [{"id": "c", "svtype": "BND", "chrom": "chr1", "pos": 100, "chr2": "chr2", "end": 1000000}]
    assert _counts(_run(tmp_path, truth, calls, size_max=50000)) == (1, 0, 0)


def test_bnd_pass_only_filter_still_applies(tmp_path):
    truth = [{"id": "t", "svtype": "BND", "chrom": "chr1", "pos": 100, "chr2": "chr2", "end": 1000000}]
    calls = [{"id": "c", "svtype": "BND", "chrom": "chr1", "pos": 100, "chr2": "chr2", "end": 1000000, "filter": "LowQual"}]
    assert _counts(_run(tmp_path, truth, calls, pass_only=True)) == (0, 0, 1)


def test_type_ignore_allows_bnd_and_tra_with_same_breakpoints(tmp_path):
    truth = [{"id": "t", "svtype": "BND", "chrom": "chr1", "pos": 1000, "chr2": "chr2", "end": 2000}]
    calls = [{"id": "c", "svtype": "TRA", "chrom": "chr1", "pos": 1000, "chr2": "chr2", "end": 2000}]
    assert _counts(_run(tmp_path, truth, calls, type_ignore=True)) == (1, 0, 0)


def test_type_ignore_does_not_mix_breakpoint_truth_with_interval_call(tmp_path):
    truth = [{"id": "t", "svtype": "TRA", "chrom": "chr1", "pos": 1000, "chr2": "chr2", "end": 1500}]
    calls = [{"id": "c", "svtype": "DEL", "chrom": "chr1", "pos": 1000, "chr2": "chr1", "end": 1500}]
    assert _counts(_run(tmp_path, truth, calls, type_ignore=True)) == (0, 1, 1)


def test_type_ignore_does_not_mix_interval_truth_with_breakpoint_call(tmp_path):
    truth = [{"id": "t", "svtype": "DEL", "chrom": "chr1", "pos": 1000, "chr2": "chr1", "end": 1500}]
    calls = [{"id": "c", "svtype": "TRA", "chrom": "chr1", "pos": 1000, "chr2": "chr2", "end": 1500}]
    assert _counts(_run(tmp_path, truth, calls, type_ignore=True)) == (0, 1, 1)
