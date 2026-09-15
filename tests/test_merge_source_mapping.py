from __future__ import annotations

from types import SimpleNamespace

from octopusv.merger.sv_merger import SVMerger
from octopusv.merger.sv_selector import select_representative_sv


FORMAT = "GT:AD:LN:ST:QV:TY:ID:SC:REF:ALT:CO"


def _make_ins(
    *,
    source_file,
    sv_id,
    pos,
    end,
    caller,
    support=10,
    qual=60,
):
    sample = {
        "GT": "0/1",
        "AD": "5,5",
        "LN": str(end - pos),
        "ST": ".",
        "QV": str(qual),
        "TY": "INS",
        "ID": sv_id,
        "SC": caller,
        "REF": "N",
        "ALT": "<INS>",
        "CO": f"chr1_{pos}-chr1_{end}",
        "original_id": sv_id,
    }

    return SimpleNamespace(
        chrom="chr1",
        pos=pos,
        sv_id=sv_id,
        ref="N",
        alt="<INS>",
        quality=str(qual),
        filter="PASS",
        info={
            "SVTYPE": "INS",
            "END": str(end),
            "SVLEN": str(end - pos),
            "SUPPORT": str(support),
        },
        format=FORMAT,
        sample=sample,
        sample_name="SAMPLE",
        source_file=str(source_file),
        start_pos=pos,
        end_pos=end,
    )


def _merge(events, input_files):
    merger = SVMerger(
        classified_events={
            "INS": {
                "chr1": events,
            }
        },
        all_input_files=input_files,
    )

    merger.merge()
    merged = merger.get_all_merged_events()

    assert len(merged) == 1

    return merger, merged[0]


def _write_one_record(
    monkeypatch,
    tmp_path,
    merger,
    event,
):
    monkeypatch.setattr(
        merger,
        "_write_vcf_header",
        lambda handle, contigs, input_files: None,
    )

    output = tmp_path / "merged.svcf"

    merger.write_results(
        output,
        [event],
        {},
        mode="caller",
        input_files=merger.all_input_files,
    )

    line = output.read_text().strip()
    fields = line.split("\t")

    info = {}

    for item in fields[7].split(";"):
        if "=" in item:
            key, value = item.split("=", 1)
            info[key] = value

    format_keys = fields[8].split(":")
    evidence_blocks = []

    for raw_sample in fields[9:]:
        values = raw_sample.split(":")

        evidence_blocks.append(
            dict(
                zip(
                    format_keys,
                    values,
                )
            )
        )

    return info, evidence_blocks


def test_exact_source_mapping_is_preserved_before_writing():
    """
    The merge layer already keeps the exact source attached to each
    evidence record. This known-good behavior must remain unchanged.
    """

    cute = "/inputs/CCS_CHM13_lra_cuteSV2_CHM13.svcf"
    debreak = "/inputs/CCS_CHM13_lra_debreak_CHM13.svcf"
    sniffles = "/inputs/CCS_CHM13_lra_sniffles2_CHM13.svcf"
    svim = "/inputs/CCS_CHM13_lra_svim_CHM13.svcf"

    group = [
        _make_ins(
            source_file=sniffles,
            sv_id="Sniffles2.INS.19BS18",
            pos=21019790,
            end=21019846,
            caller="Sniffles2_2.2",
        ),
        _make_ins(
            source_file=svim,
            sv_id="svim.INS.31319",
            pos=21019791,
            end=21019846,
            caller="SVIM-v1.4.2",
        ),
        _make_ins(
            source_file=debreak,
            sv_id="DB10911",
            pos=21019792,
            end=21019847,
            caller="DeBreak",
        ),
        _make_ins(
            source_file=cute,
            sv_id="cuteSV.INS.5830",
            pos=21019810,
            end=21019859,
            caller="cuteSV-2.0.3",
        ),
    ]

    representative = select_representative_sv(group)

    observed = [
        (
            source,
            sample_data["original_id"],
        )
        for (
            source,
            _sample_name,
            _sample_format,
            sample_data,
        ) in representative.merged_sample_records
    ]

    assert observed == [
        (
            sniffles,
            "Sniffles2.INS.19BS18",
        ),
        (
            svim,
            "svim.INS.31319",
        ),
        (
            debreak,
            "DB10911",
        ),
        (
            cute,
            "cuteSV.INS.5830",
        ),
    ]


def test_issue_189_sources_source_ids_and_evidence_are_coordered(
    monkeypatch,
    tmp_path,
):
    """
    Regression for GitHub #189.

    Input-file order deliberately differs from coordinate order.
    SOURCES[i], SOURCE_IDS[i], and evidence block i must all refer
    to the same original input record.
    """

    cute = "/inputs/CCS_CHM13_lra_cuteSV2_CHM13.svcf"
    debreak = "/inputs/CCS_CHM13_lra_debreak_CHM13.svcf"
    sniffles = "/inputs/CCS_CHM13_lra_sniffles2_CHM13.svcf"
    svim = "/inputs/CCS_CHM13_lra_svim_CHM13.svcf"

    input_files = [
        cute,
        debreak,
        sniffles,
        svim,
    ]

    events = [
        _make_ins(
            source_file=cute,
            sv_id="cuteSV.INS.5830",
            pos=21019810,
            end=21019859,
            caller="cuteSV-2.0.3",
        ),
        _make_ins(
            source_file=debreak,
            sv_id="DB10911",
            pos=21019792,
            end=21019847,
            caller="DeBreak",
        ),
        _make_ins(
            source_file=sniffles,
            sv_id="Sniffles2.INS.19BS18",
            pos=21019790,
            end=21019846,
            caller="Sniffles2_2.2",
        ),
        _make_ins(
            source_file=svim,
            sv_id="svim.INS.31319",
            pos=21019791,
            end=21019846,
            caller="SVIM-v1.4.2",
        ),
    ]

    merger, merged_event = _merge(
        events,
        input_files,
    )

    info, evidence = _write_one_record(
        monkeypatch,
        tmp_path,
        merger,
        merged_event,
    )

    assert info["SOURCES"].split(",") == [
        "CCS_CHM13_lra_cuteSV2_CHM13",
        "CCS_CHM13_lra_debreak_CHM13",
        "CCS_CHM13_lra_sniffles2_CHM13",
        "CCS_CHM13_lra_svim_CHM13",
    ]

    expected_ids = [
        "cuteSV.INS.5830",
        "DB10911",
        "Sniffles2.INS.19BS18",
        "svim.INS.31319",
    ]

    assert (
        info["SOURCE_IDS"].split(",")
        == expected_ids
    )

    assert (
        [block["ID"] for block in evidence]
        == expected_ids
    )


def test_issue_190_no_contributing_evidence_is_silently_dropped(
    monkeypatch,
    tmp_path,
):
    """
    Regression for GitHub #190.

    One source contributes two records to the same merged group.
    All five original evidence records must survive writing.
    """

    cute = "/inputs/CCS_CHM13_lra_cuteSV2_CHM13.svcf"
    debreak = "/inputs/CCS_CHM13_lra_debreak_CHM13.svcf"
    sniffles = "/inputs/CCS_CHM13_lra_sniffles2_CHM13.svcf"
    svim = "/inputs/CCS_CHM13_lra_svim_CHM13.svcf"

    input_files = [
        cute,
        debreak,
        sniffles,
        svim,
    ]

    events = [
        _make_ins(
            source_file=cute,
            sv_id="cuteSV.INS.22",
            pos=1202256,
            end=1202329,
            caller="cuteSV-2.0.3",
            support=16,
            qual=152.7,
        ),
        _make_ins(
            source_file=debreak,
            sv_id="debreak.DB56",
            pos=1202242,
            end=1202319,
            caller="DeBreak",
        ),
        _make_ins(
            source_file=sniffles,
            sv_id="Sniffles2.INS.33S0",
            pos=1202219,
            end=1202294,
            caller="Sniffles2_2.2",
        ),
        _make_ins(
            source_file=svim,
            sv_id="svim.INS.57",
            pos=1202145,
            end=1202220,
            caller="SVIM-v1.4.2",
        ),
        _make_ins(
            source_file=svim,
            sv_id="svim.INS.58",
            pos=1202237,
            end=1202276,
            caller="SVIM-v1.4.2",
        ),
    ]

    merger, merged_event = _merge(
        events,
        input_files,
    )

    _info, evidence = _write_one_record(
        monkeypatch,
        tmp_path,
        merger,
        merged_event,
    )

    observed_ids = {
        block["ID"]
        for block in evidence
    }

    assert observed_ids == {
        "cuteSV.INS.22",
        "debreak.DB56",
        "Sniffles2.INS.33S0",
        "svim.INS.57",
        "svim.INS.58",
    }