# FLORA-seq data analysis

This repository contains the analysis workflows used for the FLORA-seq
manuscript. Code is grouped by biological dataset; numbered files follow the
approximate analysis order within each dataset.

## Repository structure

```text
1-cell-line/              Cell-line mixing and cross-platform benchmarks
2-CAF/                    Breast cancer fibroblast analysis
3-mHypo/                  Mouse hypothalamus analysis
config.example.sh          Example path configuration
```

## Path configuration

The configurable scripts use three environment variables:

```bash
cp config.example.sh config.sh
# edit config.sh, then
source config.sh
```

- `FLORASEQ_DATA_ROOT`: input matrices, FASTQ/BAM files and model outputs.
- `FLORASEQ_REFERENCE_ROOT`: genome, GTF and BED reference files.
- `FLORASEQ_RESULTS_ROOT`: analysis outputs.

If these variables are not set, the configurable scripts use `data/`,
`data/reference/` and `results/` under the repository root. Individual scripts
also expose command-line arguments or analysis-specific environment variables
for overriding inputs.

Notebook-style scripts use generic `/path/to/...` placeholders for external
inputs and software. Replace these placeholders or use the corresponding
environment variables before running the workflows.

## Software

The workflows use Python/R single-cell packages together with command-line
tools including Cell Ranger, STAR/STARsolo, zUMIs, samtools, bedtools and
RSeQC. Exact requirements depend on the module; see the dataset-level README
and script headers before running.
