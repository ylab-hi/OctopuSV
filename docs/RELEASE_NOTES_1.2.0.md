# OctopuSV 1.2.0

This release makes caller-mode merging much faster and lighter, fixes how sample-mode tools treat samples that carry evidence without an ALT genotype, and lets sample mode accept callers that do not name themselves. The SVCF format remains at version 1.1.

## What's changed

### Performance

- **Caller-mode merging in partitions:** Caller-mode merging now processes one partition at a time, reusing the disk-backed spooling introduced for sample mode in 1.1.0 ([#194](https://github.com/ylab-hi/OctopuSV/issues/194)). On a 6-caller union of 315,491 merged records, run time fell from 197 s to 68 s and peak memory from 5.07 GB to 1.06 GB. The CLI is unchanged. `--upsetr` still runs in memory.
- **Faster sample-mode lookups:** Input paths and sample labels are now resolved once per run instead of once per record. A 500-sample simulated cohort merges in about 3 minutes with 0.35 GB peak memory on one CPU.

Both changes produce byte-identical output (apart from `##fileDate`) on our regression sets: a 6-caller union, a 6-sample merge and the 500-sample simulated cohort.

### Sample mode

- **Callers that do not name themselves:** Some callers (for example SVDSS) write no `##source`, so their records carry `SC=.`. A sample whose evidence in a merged event is one such record now counts as one voter for that sample instead of stopping the merge. If one sample of one merged event mixes named and unnamed evidence, OctopuSV stops with an error that names the input file and record IDs, rather than guessing which caller the unnamed evidence belongs to.

### `subset` and `inspect`

- **Evidence is no longer mistaken for padding:** A sample block was treated as an empty placeholder whenever its GT was `0/0` or `./.`, even when it held real evidence (for example SVDSS calls with GT `0/0`, or SVIM calls with GT `./.`). Sample-mode `subset` could therefore drop records whose only evidence in the selected samples had such a genotype, and `inspect` could report false warnings. Both now follow the SVCF 1.1 definition: a sample block is a placeholder only when it records no observation (`UC=0`, `UV=0`) and carries no `ID` or `SC`.
- **Carrier counts in `inspect`:** Carriers are now counted with the same genotype rules used for consensus voting, so phased (`0|0`), haploid (`0`) and missing phased (`.|.`) genotypes are classified correctly. The JSON report adds `carrier_samples` and `evidence_samples`. The SOURCES check is now reported as `W_SAMPLE_SOURCES_EVIDENCE_MISMATCH` and compares SOURCES with the samples that carry evidence.

### Representative selection

- **Non-integer SUPPORT values:** When an input `SUPPORT` value is not a non-negative integer (SVCF 1.1 allows `.` or a non-negative integer), it is treated as unknown during representative selection, with one warning per input file, instead of stopping the merge. Ties are broken on the scored SUPPORT and QUAL values.

## Validation

We tested the changes with the full regression suite, the regression sets listed above, and a read-only correctness audit: a `correct` rerun reproduced its earlier output exactly, `svcf2vcf` output checks passed for caller-mode, sample-mode and cohort merges, and in the simulated cohort, where the true origin of every record is known, no input evidence was lost and no merge or split broke the grouping rules.

We also reran the real 10-caller VCaP pipeline (`correct` → `filter` → `merge --union`). Every difference from an earlier run of the same pipeline traced back to changes documented in 1.1.0, and every filtered input record appears in the merged union. The 10-caller merge took 29 s with 0.3 GB peak memory.

## Upgrading from 1.1.0

**No rerun of `correct` is needed:** `correct` output is unchanged. In our regression tests, merge output is unchanged for inputs that 1.1.0 already accepted. If you used `subset` on sample-mode SVCFs that contain `0/0` or `./.` evidence, rerun it.

**Sample-mode `SOURCE_IDS`:** Each sample contributes one representative ID to `SOURCE_IDS`. When a sample has several records that merge into the same event (for example, near-duplicate calls from one caller), all of them contribute to that sample's genotype, but only one ID is listed.

**SVCF remains version 1.1.** Its required fields and FORMAT layouts have not changed; see [`SVCF_specifications.md`](SVCF_specifications.md).

**Known limitation:** Representative selection still uses caller QUAL values and caller-priority labels. QUAL scales are not directly comparable across callers. The scoring model is unchanged in this release.

For the previous release, see [OctopuSV 1.1.0](RELEASE_NOTES_1.1.0.md).

## Citation

Guo Q, Li Y, Wang T-Y, Ramakrishnan A, Yang R. **OctopuSV and TentacleSV: a one-stop toolkit for multi-sample, cross-platform structural variant comparison and analysis.** *Bioinformatics*. 2025. btaf599. https://doi.org/10.1093/bioinformatics/btaf599
