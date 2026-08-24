#!/usr/bin/env bash
# =============================================================================
# Module : 2-genebody / FLORA-seq genebody coverage by transcript length (3 kb)
# Dataset: FLORA-seq 3' (mix3) and Random (mixR), human-mapped BAMs
# Source : 1_Cell_line/0-bin/10-coverage-summary.sh
#
# Final figure used the 3 kb transcript-length split:
#   short: < 3000 bp   -> RseQC_floraseq_3k_short
#   long:  >= 3000 bp  -> RseQC_floraseq_3k_long
#
# FLORA-seq BAMs (listed in floraseq_bam.txt):
#   .../mix3_human/outs/possorted_genome_bam_filter_sub10.bam
#   .../mixR_human/outs/possorted_genome_bam_filter.bam
# =============================================================================

# =============================================================================
# 1. Split reference BED at 3 kb (THRESHOLDS = [3000] in 3-split_bed_by_length.py)
# =============================================================================
cd /path/to/reference/RSeQC_ref/
python 3-split_bed_by_length.py


# =============================================================================
# 2. RSeQC geneBody_coverage.py — short (< 3000 bp)
# =============================================================================
cd /path/to/floraseq/benchmark/coverage

/path/to/software/miniforge3/envs/sc-bio/bin/geneBody_coverage.py \
  -r /path/to/reference/RSeQC_ref/hg38_GENCODE.v38_lt3000.bed \
  -i floraseq_bam.txt \
  -o RseQC_floraseq_3k_short


# =============================================================================
# 3. RSeQC geneBody_coverage.py — long (>= 3000 bp)
# =============================================================================
/path/to/software/miniforge3/envs/sc-bio/bin/geneBody_coverage.py \
  -r /path/to/reference/RSeQC_ref/hg38_GENCODE.v38_ge3000.bed \
  -i floraseq_bam.txt \
  -o RseQC_floraseq_3k_long


# =============================================================================
# 4. Merge tables and plot
# =============================================================================
# Outputs:
#   genebody_coverage_floraseq_3k_len.csv
#   genebody_coverage_floraseq_3k_len.pdf
python 5-plot_coverage_by_length.py
