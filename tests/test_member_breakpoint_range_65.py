"""#6.5: record reported member breakpoints, without touching SV grouping."""

from types import SimpleNamespace

from typer.testing import CliRunner

from octopusv.cli.cli import app
from octopusv.merger.member_breakpoint_range import (
    MEMBER_POS_RANGE, MEMBER_END_RANGE, add_member_breakpoint_ranges,
)
from octopusv.merger.merge_split import GroupSplitter, select_with_split


def _event(svtype, chrom, pos, mate_chrom, mate_pos, *, extra=None):
    return SimpleNamespace(
        chrom=chrom, pos=pos, start_chrom=chrom, start_pos=pos,
        end_chrom=mate_chrom, end_pos=mate_pos, sv_type=svtype,
        info={"SVTYPE": svtype, "END": str(mate_pos), "CHR2": mate_chrom, **(extra or {})},
    )


def test_tra_flipped_members_align_to_representative_chromosome_sides():
    rep = _event("TRA", "chr17", 75792754, "chr16", 84992411)
    members = [rep,
               _event("TRA", "chr16", 84992415, "chr17", 75792731),
               _event("TRA", "chr17", 75792754, "chr16", 84992411)]
    add_member_breakpoint_ranges(rep, members)
    assert rep.info[MEMBER_POS_RANGE] == "-23,0"
    assert rep.info[MEMBER_END_RANGE] == "0,4"


def test_simple_del_two_breakpoints():
    rep = _event("DEL", "chr1", 1000, "chr1", 1300)
    other = _event("DEL", "chr1", 995, "chr1", 1310)
    add_member_breakpoint_ranges(rep, [rep, other])
    assert rep.info[MEMBER_POS_RANGE] == "-5,0"
    assert rep.info[MEMBER_END_RANGE] == "0,10"


def test_ins_has_only_one_physical_breakpoint():
    rep = _event("INS", "chr1", 1000, "chr1", 1200, extra={MEMBER_END_RANGE:"5,20"})
    other = _event("INS", "chr1", 1003, "chr1", 6000)
    add_member_breakpoint_ranges(rep, [rep, other])
    assert rep.info[MEMBER_POS_RANGE] == "0,3"
    assert MEMBER_END_RANGE not in rep.info


def test_single_member_does_not_inherit_an_old_merger_range():
    rep = _event("TRA", "chr1", 100, "chr2", 200, extra={
        MEMBER_POS_RANGE: "-100,100", MEMBER_END_RANGE: "-200,200"})
    add_member_breakpoint_ranges(rep, [rep])
    assert MEMBER_POS_RANGE not in rep.info
    assert MEMBER_END_RANGE not in rep.info


def test_sample_group_recomputes_input_rep_range_not_nested_range():
    rep = _event("DEL", "chr1", 1000, "chr1", 1300, extra={MEMBER_POS_RANGE:"-50,30"})
    other = _event("DEL", "chr1", 1007, "chr1", 1299, extra={MEMBER_POS_RANGE:"-80,100"})
    add_member_breakpoint_ranges(rep, [rep, other])
    assert rep.info[MEMBER_POS_RANGE] == "0,7"
    assert rep.info[MEMBER_END_RANGE] == "-1,0"


def test_incompatible_member_endpoints_do_not_produce_partial_range():
    rep = _event("TRA", "chr1", 100, "chr2", 200, extra={MEMBER_POS_RANGE:"-8,9"})
    bad = _event("TRA", "chr1", 102, "chr3", 201)
    add_member_breakpoint_ranges(rep, [rep, bad])
    assert MEMBER_POS_RANGE not in rep.info
    assert MEMBER_END_RANGE not in rep.info


def test_same_positions_two_real_members_can_report_zero_zero():
    rep = _event("INV", "chr3", 100, "chr3", 400)
    add_member_breakpoint_ranges(rep, [rep, _event("INV", "chr3", 100, "chr3", 400)])
    assert rep.info[MEMBER_POS_RANGE] == "0,0"
    assert rep.info[MEMBER_END_RANGE] == "0,0"


def test_select_with_split_generates_ranges_after_final_split():
    a = _event("DEL", "chr1", 100, "chr1", 155)
    b = _event("DEL", "chr1", 103, "chr1", 157)
    c = _event("DEL", "chr1", 101, "chr1", 212)
    for i, ev in enumerate((a,b,c)):
        # #192 only splits non-equivalent structures from the SAME source.
        ev.source_file = "/caller_a.vcf" if i != 1 else "/caller_b.vcf"
        ev.sv_id = f"call{i}"
        ev.info['SVLEN'] = str(ev.end_pos - ev.pos)
    representatives = select_with_split([[a,b,c]], GroupSplitter("caller"), lambda part:part[0])
    assert len(representatives)==2
    assert representatives[0].info[MEMBER_POS_RANGE]=="0,3"
    assert representatives[0].info[MEMBER_END_RANGE]=="0,2"
    assert MEMBER_POS_RANGE not in representatives[1].info
    assert representatives[0].info['MERGE_SPLIT']==representatives[1].info['MERGE_SPLIT']


def test_vcap_sv3_range_is_written_and_survives_standard_vcf_export(tmp_path):
    from pathlib import Path
    fixture = Path(__file__).parent / 'data/merge_guard_vcap/tra_sv3_zdhhc7_8_callers/inputs'
    src = fixture / 'manifest.tsv'
    merged = tmp_path / 'merged.svcf'
    runner = CliRunner()
    result = runner.invoke(app, ['merge', '--mode', 'caller', '--union',
                                 '--input-list', str(src), '-o', str(merged)])
    assert result.exit_code == 0, result.output
    text = merged.read_text()
    assert text.count('##INFO=<ID=MERGE_POS_RANGE,') == 1
    assert text.count('##INFO=<ID=MERGE_END_RANGE,') == 1
    lines = [line for line in text.splitlines() if '\tBND00022166\t' in line]
    assert len(lines) == 1
    assert 'MERGE_POS_RANGE=-23,0' in lines[0]
    assert 'MERGE_END_RANGE=0,4' in lines[0]
    result = runner.invoke(app, ['validate-svcf', '-i', str(merged)])
    assert result.exit_code == 0, result.output
    exported = tmp_path / 'export.vcf'
    result = runner.invoke(app, ['svcf2vcf', '-i', str(merged), '-o', str(exported)])
    assert result.exit_code == 0, result.output
    vcf = exported.read_text()
    assert '##INFO=<ID=MERGE_POS_RANGE,Number=2,Type=Integer,' in vcf
    assert 'MERGE_POS_RANGE=-23,0' in vcf
    assert 'MERGE_END_RANGE=0,4' in vcf
