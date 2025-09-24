# OctopuSV: Advanced structural variant analysis toolkit 🐙

<p align="center">
  <img src="https://github.com/ylab-hi/octopusV/blob/main/imgs/logo.png" width="40%" height="40%">
</p>

[![PyPI version](https://badge.fury.io/py/octopusv.svg)](https://badge.fury.io/py/octopusv)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

**OctopuSV** solves three critical problems in structural variant (SV) analysis:

1. **Smart BND standardization** - Converts paired BND records into standard SV types (DEL/INV/DUP/TRA), while preserving potential complex rearrangements as BNDs
2. **Multi-caller integration** - Merge SVs from different tools (Manta, Sniffles, PBSV, etc.) with flexible logic 
3. **Somatic variant calling** - Extract tumor-specific SVs by comparing tumor vs normal samples

Whether you're analyzing single samples, cohorts, or tumor/normal pairs, OctopuSV standardizes your workflow from raw calls to publication-ready results.

## Key Features

* **One-command somatic calling**: Extract tumor-specific SVs with `octopusv somatic`
* **Flexible SV merging**: Boolean logic, intersection, union, or custom expressions across samples/callers
* **BND standardization**: Converts complex breakend notation to standard SV types
* **Cross-platform support**: Works with Illumina, PacBio, ONT data from 10+ popular callers
* **Built-in benchmarking**: Compare against truth sets with precision/recall metrics
* **Rich visualizations**: Interactive HTML reports and publication-ready plots

## Installation

```bash
pip install octopusv
```

---

## Quick Start

### 1. Correct Ambiguous BND Annotations

```bash
# Basic correction
octopusv correct input.vcf output.vcf

# With position tolerance control
octopusv correct -i input.vcf -o output.vcf --pos-tolerance 5

# Apply quality filters
octopusv correct -i input.vcf -o output.vcf --min-svlen 50 --max-svlen 100000 --filter-pass
```

### 2. Merge SV Calls (Multi-caller or Multi-sample)

```bash
# Basic intersection: SVs found by ALL callers
octopusv merge -i manta.svcf sniffles.svcf pbsv.svcf -o intersection.svcf --intersect

# Union: SVs found by ANY caller  
octopusv merge -i caller1.svcf caller2.svcf caller3.svcf -o union.svcf --union

# Specific caller: SVs unique to one caller
octopusv merge -i manta.svcf sniffles.svcf -o manta_specific.svcf --specific manta.svcf

# Minimum support: SVs supported by at least N callers
octopusv merge -i a.svcf b.svcf c.svcf d.svcf -o supported.svcf --min-support 3

# Complex Boolean logic: (A AND B) but NOT (C OR D)
octopusv merge -i A.svcf B.svcf C.svcf D.svcf \
  --expression "(A AND B) AND NOT (C OR D)" -o filtered.svcf

# Multi-sample mode with custom names
octopusv merge -i sample1.svcf sample2.svcf sample3.svcf \
  --mode sample --sample-names Patient1,Patient2,Patient3 \
  --min-support 2 -o cohort.svcf

# Generate intersection plot
octopusv merge -i a.svcf b.svcf c.svcf -o merged.svcf --intersect \
  --upsetr --upsetr-output venn_diagram.png
```

<p align="center">
  <img src="https://github.com/ylab-hi/octopusV/blob/main/imgs/up_upset.png" width="70%" height="70%">
</p>

### 3. Somatic SV Calling (NEW!)

Extract tumor-specific structural variants by comparing tumor and normal samples:
```bash
# Basic somatic calling
octopusv somatic -t tumor.svcf -n normal.svcf -o somatic.svcf

# With custom matching parameters
octopusv somatic -t tumor.svcf -n normal.svcf -o somatic.svcf \
  --max-distance 100 --min-jaccard 0.8

# Convert to standard VCF for downstream analysis
octopusv svcf2vcf -i somatic.svcf -o somatic.vcf
```

### 4. Benchmark Against Truth Sets

```bash
octopusv benchmark truth.vcf calls.svcf \
  -o benchmark_results \
  --reference-distance 500 \
  --size-similarity 0.7 \
  --reciprocal-overlap 0.0 \
  --size-min 50 --size-max 50000
```

### 5. Generate Statistics and Visualizations

```bash
# Basic stat collection
octopusv stat -i input.svcf -o stats.txt

# Add HTML report
octopusv stat -i input.svcf -o stats.txt --report

# Plot figures from stats
octopusv plot stats.txt -o figure_prefix
```

The `--report` flag outputs an interactive HTML report:

* SV type and size distributions
* Chromosome breakdowns
* Quality score summaries
* Genotype and depth features

<p align="center">
  <img src="https://github.com/ylab-hi/octopusV/blob/main/imgs/html_example.png" width="70%" height="70%">
</p>

### 6. Format Conversion

```bash
# To BED
octopusv svcf2bed -i input.svcf -o output.bed

# To BEDPE
octopusv svcf2bedpe -i input.svcf -o output.bedpe

# To standard VCF
octopusv svcf2vcf -i input.svcf -o output.vcf
```

---

## Example Visualizations

OctopusV generates publication-ready visualizations:

### Chromosome Distribution

<p align="center">
  <img src="https://github.com/ylab-hi/octopusV/blob/main/imgs/chromosome_distribution.png" width="50%" height="50%">
</p>

### SV Type Distribution

<p align="center">
  <img src="https://github.com/ylab-hi/octopusV/blob/main/imgs/sv_types.png" width="50%" height="50%">
</p>

### SV Size Distribution

<p align="center">
  <img src="https://github.com/ylab-hi/octopusV/blob/main/imgs/sv_sizes.png" width="50%" height="50%">
</p>

---

See the companion pipeline: [TentacleSV](https://github.com/ylab-hi/TentacleSV)

---

## 🧪 Citation

If you use **OctopuSV**, please cite:

> Guo Q, Li Y, Wang T, Ramakrishnan A, Yang R. *OctopuSV and TentacleSV: a one-stop toolkit for multi-sample, cross-platform structural variant comparison and analysis*. bioRxiv. 2025. doi: [10.1101/2025.03.24.645012](https://doi.org/10.1101/2025.03.24.645012)

```bibtex
@article{guo2025octopusv,
  title={OctopuSV and TentacleSV: a one-stop toolkit for multi-sample, cross-platform structural variant comparison and analysis},
  author={Guo, Qingxiang and Li, Yangyang and Wang, Tingyou and Ramakrishnan, Abhi and Yang, Rendong},
  journal={bioRxiv},
  year={2025},
  publisher={Cold Spring Harbor Laboratory},
  doi={10.1101/2025.03.24.645012},
  url={https://www.biorxiv.org/content/10.1101/2025.03.24.645012v1}
}
```

---

## Contributing

We welcome issues, suggestions, and pull requests!

```bash
git clone https://github.com/ylab-hi/OctopuSV.git
cd OctopuSV
poetry install
pre-commit run -a
```

## Contact

* GitHub Issues: [https://github.com/ylab-hi/octopusV/issues](https://github.com/ylab-hi/octopusV/issues)
* Email: [qingxiang.guo@northwestern.edu](mailto:qingxiang.guo@northwestern.edu)
* Email: [yangyang.li@northwestern.edu](mailto:yangyang.li@northwestern.edu)
