# OctopuSV 1.1.0

This release fixes several issues found in real SV callsets and makes multi-sample merging less memory-intensive. The SVCF format remains at version 1.1.

## What's changed

### VCF parsing and normalization

- **Caller-specific support counts:** `correct` now reads support from the appropriate fields for each caller. A reported zero is kept distinct from a missing value.
- **BND/TRA handling:** Improved mate validation and conversion preserve distinct breakend junctions. Reciprocal `MATEID` records are checked against breakpoint coordinates and orientation.
- **Length filtering ([#193](https://github.com/ylab-hi/OctopuSV/issues/193)):** `--min-svlen` no longer removes TRA/BND calls simply because a caller reports `SVLEN=0`.
- **Caller compatibility:** Fixed SvABA sample-column parsing and recovered DELLY caller names from `SVMETHOD` when `##source` is absent. DeBreak TRA records without usable remote coordinates are now skipped with a separate report.

### Merging and output

- **Distinct SV structures ([#192](https://github.com/ylab-hi/OctopuSV/issues/192)):** Different DEL/DUP/INS structures from the same source are no longer collapsed into one merged event. Related split records are marked with `MERGE_SPLIT`.
- **Lower memory use ([#194](https://github.com/ylab-hi/OctopuSV/issues/194)):** Sample-mode merging now processes SVs in disk-backed partitions. The CLI is unchanged, and regression tests confirm the same output as the previous in-memory path.
- **Member breakpoint ranges:** Optional `MERGE_POS_RANGE` and `MERGE_END_RANGE` INFO fields summarize the coordinates reported by members of a merged event. These are observed coordinate ranges, **not confidence intervals**.
- **TRA option ([#191](https://github.com/ylab-hi/OctopuSV/issues/191)):** `--no-tra-strand-consistency` now works as intended.

## Validation

We tested the changes with the full regression suite and real 10-caller VCaP data. The member-range update left the existing merged output unchanged apart from its new INFO annotations.

In a sample-mode stress test using **250 simulated samples built from real caller SVCFs**, OctopuSV processed **6.1 million input records** in about **29 minutes**, with **7.1 GiB peak memory** on one CPU. The output passed SVCF validation. This was a simulated cohort, not 250 independent biological samples.

## Upgrading from 1.0.0

We recommend rerunning **`correct` → `merge` from the original VCFs** to pick up the normalization and caller-identity fixes. This is particularly important for DELLY files that previously had `SC=.`.

**SVCF remains version 1.1.** Its required fields and FORMAT layouts have not changed. The new merge annotations are optional; see [`SVCF_specifications.md`](SVCF_specifications.md).

**Known limitation:** Representative selection still uses caller QUAL values and caller-priority labels. QUAL scales are not directly comparable across callers, and some program labels are not recognized by the current priority mapping. The scoring model is unchanged in this release.

For the previous release, see [OctopuSV 1.0.0](RELEASE_NOTES_1.0.0.md).

## Citation

Guo Q, Li Y, Wang T-Y, Ramakrishnan A, Yang R. **OctopuSV and TentacleSV: a one-stop toolkit for multi-sample, cross-platform structural variant comparison and analysis.** *Bioinformatics*. 2025. btaf599. https://doi.org/10.1093/bioinformatics/btaf599
