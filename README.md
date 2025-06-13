# OctopuSV: Advanced structural variant analysis toolkit 🐙

<p align="center">
  <img src="https://github.com/ylab-hi/octopusV/blob/main/imgs/logo.png" width="40%" height="40%">
</p>

[![PyPI version](https://badge.fury.io/py/octopusv.svg)](https://badge.fury.io/py/octopusv)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

OctopuSV is a comprehensive toolkit for analyzing, comparing, and integrating structural variant (SV) calls from multiple platforms and callers. It addresses key challenges in SV analysis by standardizing ambiguous breakend (BND) annotations, providing flexible Boolean-based merging operations, and offering an integrated framework for SV benchmarking and visualization.

## Key Features

- **BND Correction**: Accurately convert ambiguous breakend (BND) notations into canonical SV types (INV, DUP, TRA)
- **Flexible Merging**: Combine SV calls using advanced Boolean operations, including difference and complement sets
- **Multi-caller Integration**: Integrate SV calls from various platforms (Illumina, PacBio, ONT) and callers (Manta, LUMPY, SvABA, DELLY, PBSV, Sniffles, etc.)
- **Benchmarking**: Compare SV calls against truth sets using GIAB standards with type-specific metrics
- **Visualization**: Generate publication-ready plots for comprehensive SV analysis

## TentacleSV: Automated SV Analysis Pipeline 👾 

OctopusV is paired with **TentacleSV**, a Snakemake-based workflow that automates the entire SV analysis process from raw reads to final merged callsets. TentacleSV coordinates read mapping, multi-caller variant detection, and OctopusV's correction and merging operations through a simple configuration file.

Visit [TentacleSV repository](https://github.com/ylab-hi/TentacleSV) for more information.

## Installation

```bash
pip install octopusv
```

## Quick Start

### 1. Correct Ambiguous BND Annotations

Convert breakend (BND) annotations from various SV callers to canonical SV types:

```bash
# Basic correction
octopusv correct input.vcf output.vcf

# With adjusted position tolerance
octopusv correct -i input.vcf -o output.vcf --pos-tolerance 5
```

### 2. Merge SV Calls

Combine SV calls with various strategies:

```bash
# Intersect calls from multiple callers
octopusv merge -i caller1.vcf caller2.vcf -o merged.vcf --intersect

# Union of calls
octopusv merge -i caller1.vcf caller2.vcf -o merged.vcf --union

# Get SVs supported by at least 2 callers
octopusv merge -i caller1.vcf caller2.vcf caller3.vcf -o merged.vcf --min-support 2

# Extract caller-specific SVs
octopusv merge -i caller1.vcf caller2.vcf -o specific.vcf --specific caller1.vcf

# Complex merging with Boolean expression
octopusv merge -i caller1.vcf caller2.vcf caller3.vcf -o merged.vcf \
    --expression "(caller1 AND caller2) AND NOT caller3"

# Generate UpSet plot of intersections
octopusv merge -i caller1.vcf caller2.vcf caller3.vcf -o merged.vcf \
    --intersect --upsetr --upsetr-output intersections.png
```

<p align="center">
  <img src="https://github.com/ylab-hi/octopusV/blob/main/imgs/up_upset.png" width="70%" height="70%">
</p>

### 3. Benchmark Against Truth Sets

Evaluate SV calls against a truth set using GIAB standards:

```bash
octopusv benchmark \
    truth_set.vcf \
    test_calls.vcf \
    -o benchmark_results \
    --reference-distance 500 \
    --size-similarity 0.7 \
    --reciprocal-overlap 0.0 \
    --size-min 50 \
    --size-max 50000
```

### 4. Generate Statistics and Visualizations

Analyze SV characteristics with comprehensive statistics:

```bash
# Basic statistics output
octopusv stat -i variants.vcf -o stats.txt

# With size filters
octopusv stat -i variants.vcf -o stats.txt --min-size 50 --max-size 10000

# Generate comprehensive HTML report
octopusv stat -i variants.vcf -o stats.txt --report

# Generate plots from statistics
octopusv plot stats.txt -o plots_prefix
```

The `--report` option generates an interactive HTML report that includes:
- Summary statistics for all SV types
- SV type distribution visualization
- Size distribution analysis
- Chromosome coverage plots
- Quality metrics visualization

<p align="center">
  <img src="https://github.com/ylab-hi/octopusV/blob/main/imgs/html_example.png" width="70%" height="70%">
</p>

### 5. Format Conversion

Convert between different SV representation formats:

```bash
# Convert to BED format
octopusv svcf2bed -i variants.vcf -o variants.bed

# Convert to BEDPE format
octopusv svcf2bedpe -i variants.vcf -o variants.bedpe

# Convert to standard VCF format
octopusv svcf2vcf -i variants.vcf -o variants.vcf
```

## Advanced Usage

### Merging Configuration Options

- `--max-distance`: Maximum distance between start/end positions (default: 50)
- `--max-length-ratio`: Maximum ratio between event lengths (default: 1.3)
- `--min-jaccard`: Minimum Jaccard index for overlap (default: 0.7)
- `--tra-delta`: Position uncertainty for translocations (default: 50)
- `--tra-min-overlap`: Minimum overlap ratio for translocations (default: 0.5)
- `--tra-strand-consistency`: Require strand consistency (default: True)

### Benchmark Configuration Options

- `--reference-distance`: Maximum reference location distance (default: 500bp)
- `--size-similarity`: Minimum size similarity ratio (default: 0.7)
- `--reciprocal-overlap`: Minimum reciprocal overlap (default: 0.0)
- `--type-ignore`: Allow type mismatches
- `--enable-sequence-comparison`: Enable sequence similarity comparison
- `--sequence-similarity`: Minimum sequence similarity threshold (default: 0.7)

## Example Visualizations

OctopusV creates publication-ready visualizations for comprehensive SV analysis:

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

## Research Applications

OctopusV has been designed to address critical challenges in SV analysis, particularly for:

- Standardizing ambiguous breakend (BND) annotations from diverse SV callers
- Enabling precise identification of cohort-specific SVs without custom scripting
- Automating complex SV analysis workflows for large-scale studies
- Improving precision and recall in multi-platform, multi-caller SV analysis

For more details, see our publication: [OctopusV and TentacleSV: a one-stop toolkit for multi-sample, cross-platform structural variant comparison and analysis](https://github.com/ylab-hi/OctopuSV).

## Contributing

We welcome contributions! Here's how you can help:

1. Clone the repository:

   ```bash
   git clone https://github.com/ylab-hi/OctopuSV.git
   ```

2. Create a feature branch:

   ```bash
   git checkout -b feature/your-feature-name
   ```

3. Install development dependencies:

   ```bash
   poetry install
   ```

4. Make your changes and ensure all checks pass:

   ```bash
   pre-commit run -a
   ```

5. Submit a pull request with a clear description of your changes

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Contact

For questions and feedback:

- GitHub Issues: [https://github.com/ylab-hi/octopusV/issues](https://github.com/ylab-hi/octopusV/issues)
- Email: qingxiang.guo@northwestern.edu
- Email: yangyang.li@northwestern.edu


## 📄 Citation

If you use **OctopuSV** or **TentacleSV**, please cite:

> Guo Q, Li Y, Wang T, Ramakrishnan A, Yang R. *OctopusV and TentacleSV: a one-stop toolkit for multi-sample, cross-platform structural variant comparison and analysis*. bioRxiv. 2025. doi: [10.1101/2025.03.24.645012](https://doi.org/10.1101/2025.03.24.645012)

## 📚 BibTeX

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
