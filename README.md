# OctopuSV: Advanced structural variant analysis toolkit 🐙

<p align="center">
  <img src="https://github.com/ylab-hi/octopusV/blob/main/imgs/logo.png" width="40%" height="40%">
</p>

[![PyPI version](https://badge.fury.io/py/octopusv.svg)](https://badge.fury.io/py/octopusv)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

**OctopuSV** is a high-performance structural variant (SV) analysis toolkit designed to standardize ambiguous SV annotations (e.g., BNDs), flexibly integrate multiple callers across samples or platforms, and benchmark results against trusted truth sets. With support for both **single-sample** and **multi-sample** workflows, OctopuSV enables robust and scalable SV comparison, correction, and visualization in real or simulated genomic datasets.

## Key Features

* **BND Correction**: Converts ambiguous breakend (BND) records into canonical SV types (DEL, INV, DUP, TRA), with translocation subtype classification
* **Flexible Multi-sample Merging**: Boolean logic-based merge of SVs across multiple samples or callers
* **Multi-caller & Multi-platform Integration**: Works seamlessly across Illumina, PacBio, ONT callers like Manta, LUMPY, SvABA, DELLY, PBSV, Sniffles, etc.
* **Benchmarking**: Compare SVs to truth sets with precision/recall/F1 metrics using GIAB-style evaluation
* **Statistical Summaries**: Profile SV distribution, quality, and size
* **Publication-ready Visualizations**: Output interactive HTML reports and static plots

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
# Merge across callers from same sample
octopusv merge -i manta.svcf lumpy.svcf -o merged.svcf --mode caller --caller-names Manta,LUMPY --intersect

# Merge across samples
octopusv merge -i sample1.svcf sample2.svcf sample3.svcf \
  --mode sample --sample-names HG001,HG002,HG003 \
  --min-support 2 -o shared.svcf

# Complex logic: A AND B but not C
octopusv merge -i A.svcf B.svcf C.svcf \
  --expression "(A AND B) AND NOT C" -o result.svcf

# Generate UpSet plot
octopusv merge -i a.svcf b.svcf c.svcf -o merged.svcf --intersect --upsetr --upsetr-output intersection.png
```

<p align="center">
  <img src="https://github.com/ylab-hi/octopusV/blob/main/imgs/up_upset.png" width="70%" height="70%">
</p>

### 3. Benchmark Against Truth Sets

```bash
octopusv benchmark truth.vcf calls.svcf \
  -o benchmark_results \
  --reference-distance 500 \
  --size-similarity 0.7 \
  --reciprocal-overlap 0.0 \
  --size-min 50 --size-max 50000
```

### 4. Generate Statistics and Visualizations

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

### 5. Format Conversion

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

## Application Scenarios

OctopuSV was developed to address several practical needs in SV research:

* Standardizing SVs with ambiguous BND notations
* Enabling precise cohort-level comparisons (multi-sample mode)
* Supporting accurate benchmarking with real/simulated truth sets
* Integrating and comparing SVs across platforms (e.g., Illumina + ONT)
* Automating large-scale SV analysis workflows (via TentacleSV)

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
