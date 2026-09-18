# OctopuSV 1.0.0

Released: 2026-09-18

OctopuSV 1.0.0 is the first release built around the versioned SVCF 1.1 format and the full multi-caller / multi-sample workflow. This release focuses on making merged SV records easier to interpret, safer to reuse downstream, and more reliable across larger cohorts.

## Highlights

### SVCF 1.1

OctopuSV now writes versioned SVCF 1.1 files with an explicit data model:

- `caller` mode keeps source-level evidence for merged events.
- `multi` mode stores a fixed biological-sample matrix with synthesized sample calls.

For caller-mode merges, source/evidence relationships are explicit and positional:

```text
SOURCES[i] <-> SOURCE_IDS[i] <-> evidence block i
```

Multiple evidence records from the same caller are preserved. Source labels are treated as exact, case-sensitive identities.

See `docs/SVCF_specifications.md` for the full format contract.

### Multi-sample and cohort workflows

Sample-mode merging has been reworked so that biological sample columns remain stable through merge, filtering, subsetting, and VCF export.

Sample genotypes are synthesized from unique caller votes rather than raw evidence-record counts. Unobserved sample/event combinations remain distinguishable from evidence-backed `0/0` calls.

By default, `svcf2vcf` exports true unobserved samples as `./.`. Use:

```bash
--unobserved-sample-gt ref
```

to export those placeholders as `0/0` when that interpretation is appropriate for the downstream analysis.

The current workflow has also been tested on a 300-sample cohort.

### Source/evidence mapping

Several edge cases in merged source tracking were fixed and regression-tested.

In particular:

- `SOURCES`, `SOURCE_IDS`, and evidence blocks are now generated from the same ordered record list.
- duplicate source labels are preserved when one caller contributes multiple records to the same event.
- source support is based on unique sources, not the number of evidence records.
- source IDs containing colons, including Manta event IDs, are preserved without truncation.

If you rely on `SOURCES` / `SOURCE_IDS` from older merged SVCF files, we recommend rerunning the merge with 1.0.0 rather than editing old files manually.

### BND, TRA, and GRIDSS handling

BND/TRA handling has been tightened without changing the validated core merge behavior.

- paired BND records can be resolved into standard SV types when the breakpoint evidence supports it.
- interchromosomal BND records with an explicit remote mate coordinate can be represented as `TRA`.
- unresolved BND records remain valid BND records rather than being forced into another SV type.
- `TRA` can use breakend ALT notation when orientation is known, or symbolic `<TRA>` with `CHR2` and `END` when orientation is unknown.
- true GRIDSS single-breakends are skipped by default with an explicit note and count in the output header.
- `--strict-single-breakends` restores hard-fail behavior when desired.

OctopuSV does not infer a missing remote breakpoint from GRIDSS `BEALN`.

### Merge controls

`--max-distance`, `--max-length-ratio`, and `--min-jaccard` now act as real merge controls.

`--min-jaccard` applies to DEL, DUP, and INV records and remains disabled by default (`0`).

Conflicting merge strategies now fail early instead of silently ignoring one option. For example:

```bash
--intersect --min-support 3
```

is rejected because those two options describe different support rules.

### Safer SVCF to VCF conversion

`svcf2vcf` now preserves record-level event coordinates rather than borrowing endpoints from another caller's evidence block.

For insertions:

```text
SVCF: END = POS + SVLEN
VCF:  END = POS
```

The conversion also preserves multi-sample columns and keeps `SOURCES` / `SOURCE_IDS` visible in the VCF INFO field where applicable.

### Validation and metadata safety

The SVCF validator now checks the versioned SVCF 1.1 structure, including:

- mode-specific FORMAT schemas
- source/evidence positional consistency
- sample-column counts
- BND ALT / `CHR2` / `END` agreement
- TRA breakpoint representation
- required INFO fields
- reserved characters in source labels
- structural coordinate rules

Reference metadata is handled conservatively:

- different reference/assembly strings alone do not force a merge failure.
- exact input metadata can be retained with a warning when compatibility cannot be proven from the strings alone.
- known contig-length conflicts fail before merge.

`normalize-contigs` only changes naming style, such as `1` to `chr1`. It does not perform liftover.

### Compatibility updates

This release also includes compatibility and parsing improvements for:

- GRIDSS
- SVIM-ASM
- DRAGEN CNV
- multi-sample raw VCF input
- colon-containing record IDs and BND ALT strings
- custom FILTER definitions carried through merged output and VCF export

DRAGEN `SVTYPE=CNV` records are normalized to supported `DEL` / `DUP` events when the call provides enough information to determine the direction.

## Notes for upgrading from 0.4.x

SVCF 1.1 is a versioned format contract. Legacy unversioned SVCF files remain readable through compatibility paths where supported, but new output should use the SVCF 1.1 model.

Caller-mode SVCF is an intermediate format and is not guaranteed to behave like a conventional VCF sample matrix. Before using merged SVCF with tools such as bcftools or vcftools, convert it first:

```bash
octopusv svcf2vcf -i merged.svcf -o merged.vcf
```

For cohort analysis, the recommended workflow is:

```text
raw VCF per biological sample
    -> octopusv correct
    -> caller-mode SVCF per sample
    -> octopusv merge --mode sample
    -> multi-mode SVCF
```

Historical multi-sample output produced directly by `octopusv correct` is retained as a legacy compatibility format and should not be used as a substitute for the SVCF 1.1 sample-mode workflow.

## Validation

Before release, OctopuSV 1.0.0 was checked with:

- the full regression test suite
- real heterogeneous caller workflows
- caller-mode and sample-mode downstream conversion
- GRIDSS BND/single-breakend handling
- GRCh37 / GRCh38 contig-length compatibility checks
- a 300-sample cohort smoke test
- matched-input comparison against the previous production merge workflow
- a historical GitHub issue audit

The core event-grouping behavior was kept stable unless a reproducible bug required a change.

## Citation

If you use OctopuSV in your research, please cite:

Guo Q, Li Y, Wang T-Y, Ramakrishnan A, Yang R. **OctopuSV and TentacleSV: a one-stop toolkit for multi-sample, cross-platform structural variant comparison and analysis.** *Bioinformatics*. 2025. btaf599. https://doi.org/10.1093/bioinformatics/btaf599
