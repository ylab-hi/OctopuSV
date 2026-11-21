"""Sample-mode lookups done once per run, with unchanged results.

* NameMapper builds its path -> label table once (it used list.index per
  sample column of every record).
* Partitioned sample merge resolves each input path once per run (it called
  realpath for every evidence record in every bucket).
"""

import os

import pytest
from typer.testing import CliRunner

from octopusv.cli.cli import app
from octopusv.merger.name_mapper import NameMapper


def test_name_mapper_labels_order_and_first_duplicate_wins():
    files = ["/d/a.svcf", "/d/b.svcf", "/d/a.svcf", "/d/c.svcf"]
    mapper = NameMapper(files, mode="sample", custom_names=["A", "B", "A2", "C"])

    # Historical list.index behaviour: a repeated path always resolves to its
    # first occurrence.  (merge itself rejects duplicate inputs earlier.)
    assert mapper.get_all_display_names() == ["A", "B", "A", "C"]
    assert mapper.get_display_name("/d/a.svcf") == "A"
    assert mapper.get_display_name("/d/c.svcf") == "C"
    assert mapper.convert_source_string("/d/c.svcf,/d/b.svcf") == "C,B"


def test_name_mapper_unknown_path_still_refused():
    mapper = NameMapper(["/d/a.svcf"], mode="sample", custom_names=["A"])
    with pytest.raises(ValueError, match="is not present in the NameMapper input list"):
        mapper.get_display_name("/d/other.svcf")


def test_name_mapper_default_labels_and_returned_list_is_a_copy():
    mapper = NameMapper(["/d/s1.svcf.gz", "/d/s2.svcf"], mode="sample")
    names = mapper.get_all_display_names()
    assert names == ["s1", "s2"]
    names.append("changed")
    assert mapper.get_all_display_names() == ["s1", "s2"]


def _svcf(path, sample, n_records):
    lines = [
        "##fileformat=SVCFv1.1",
        "##contig=<ID=chr1,length=10000000>",
        f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample}",
    ]
    for i in range(n_records):
        pos = 1000 + i * 5000
        end = pos + 500
        rid = f"{sample}.{i}"
        lines.append(
            f"chr1\t{pos}\t{rid}\tN\t<DEL>\t30\tPASS\t"
            f"SVTYPE=DEL;END={end};SVLEN=500;SUPPORT=10\t"
            "GT:AD:LN:ST:QV:TY:ID:SC:REF:ALT:CO\t"
            f"0/1:5,10:500:+-:30:DEL:{rid}:sniffles:N:<DEL>:chr1_{pos}-chr1_{end}"
        )
    path.write_text("\n".join(lines) + "\n")


def test_partitioned_sample_merge_resolves_paths_once_per_run(tmp_path, monkeypatch):
    inputs = []
    for index in range(3):
        path = tmp_path / f"s{index}.svcf"
        _svcf(path, f"s{index}", 300)
        inputs += ["-i", str(path)]

    calls = []
    real = os.path.realpath

    def counted(*args, **kwargs):
        calls.append(args[0])
        return real(*args, **kwargs)

    monkeypatch.setattr(os.path, "realpath", counted)
    output = tmp_path / "merged.svcf"
    result = CliRunner().invoke(app, ["merge", *inputs, "-o", str(output), "--mode", "sample", "--union"])
    assert result.exit_code == 0, result.output

    rows = [line for line in output.read_text().splitlines() if not line.startswith("#")]
    assert len(rows) == 300
    # 900 evidence records: realpath is no longer called per record.
    assert len(calls) < 60, len(calls)
