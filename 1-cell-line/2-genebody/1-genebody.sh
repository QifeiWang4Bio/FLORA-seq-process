#!/usr/bin/env bash
# =============================================================================
# Module : 2-genebody / Cross-technology gene body coverage (RSeQC)
# Dataset: Cell-line benchmark (human-mapped BAMs, hg38)
# Source : 1_Cell_line/0-bin/3-genebody.sh
#
# Workflow:
#   1) Filter unmapped / secondary / supplementary reads
#   2) Subsample BAMs to comparable read depth across technologies
#   3) Run RSeQC geneBody_coverage.py (see SLURM examples below)
#   4) Merge tables and plot (2-plot_genebody.py)
#
# NOTE
# - Absolute HPC paths from the original analysis; edit for your environment.
# - Sections are meant to be run selectively, not end-to-end.
#
# Environment (typical):
#   conda activate sc-bio
#   samtools, RSeQC (geneBody_coverage.py), subset-bam (optional)
# =============================================================================

# =============================================================================
# Design note (BAM inventory used for genebody coverage)
# =============================================================================
# For genebody coverage, we did not restrict to a single cell type. Instead, we
# used the BAM files produced by each technology (all human cells mapped to
# hg38), randomly subsampled each BAM so that read counts were kept at a
# comparable scale across technologies, and then computed genebody coverage.


# =============================================================================
# 0. Optional: subset HEK293T cells (not required for genebody design above)
# =============================================================================
# subset-bam \
#   --bam possorted_genome_bam.bam \
#   --cell-barcodes HEK293T_barcodes.txt \
#   --out-bam HEK293T_only.bam


# =============================================================================
# 1. Filter BAMs and subsample to comparable depth
#    Flags: -F 4 (unmapped), -F 256 (secondary), -F 2048 (supplementary)
# =============================================================================

# --- 1.1 FLORA-seq 3' only (mix3) ---
cd /path/to/floraseq/cell_line/mix3_human/outs/
samtools view -b -F 4 -F 256 -F 2048 \
  possorted_genome_bam.bam > possorted_genome_bam_filter.bam
samtools index possorted_genome_bam_filter.bam

# subsample 10%
samtools view -b -s 0.1 possorted_genome_bam_filter.bam > possorted_genome_bam_filter_sub10.bam
samtools index possorted_genome_bam_filter_sub10.bam

# --- 1.2 FLORA-seq Random only (mixR); used as baseline (no subsample) ---
cd /path/to/floraseq/cell_line/mixR_human/outs/
samtools view -b -F 4 -F 256 -F 2048 \
  possorted_genome_bam.bam > possorted_genome_bam_filter.bam
samtools index possorted_genome_bam_filter.bam

# --- 1.3 FLORA-seq 3'+Random (mixR3) ---
cd /path/to/floraseq/cell_line/mixR3_human/outs/
samtools view -b -F 4 -F 256 -F 2048 \
  possorted_genome_bam.bam > possorted_genome_bam_filter.bam
samtools index possorted_genome_bam_filter.bam

# subsample 7%
samtools view -b -s 0.07 possorted_genome_bam_filter.bam > possorted_genome_bam_filter_sub7.bam
samtools index possorted_genome_bam_filter_sub7.bam

# --- 1.4 VASA-seq ---
cd /path/to/floraseq/benchmark/vasaseq/star
samtools view -b -F 4 -F 256 -F 2048 \
  vasa_293tAligned.sortedByCoord.out.bam > vasa_293tAligned.sortedByCoord.out_filter.bam
samtools index vasa_293tAligned.sortedByCoord.out_filter.bam

# subsample 5%
samtools view -b -s 0.05 vasa_293tAligned.sortedByCoord.out_filter.bam \
  > vasa_293tAligned.sortedByCoord.out_sub5.bam
samtools index vasa_293tAligned.sortedByCoord.out_sub5.bam

# --- 1.5 10x 3' ---
cd /path/to/floraseq/benchmark/10x_3/293t_mix_1k_v3-1/293t_mix_1k_v3-1_human/outs
samtools view -b -F 4 -F 256 -F 2048 \
  possorted_genome_bam.bam > possorted_genome_bam_filter.bam
samtools index possorted_genome_bam_filter.bam

# subsample 17%
samtools view -b -s 0.17 possorted_genome_bam_filter.bam > possorted_genome_bam_filter_sub17.bam
samtools index possorted_genome_bam_filter_sub17.bam

# --- 1.6 10x 5' ---
cd /path/to/floraseq/benchmark/10x_5/10k_hgmm_5pv2_human/outs
samtools view -b -F 4 -F 256 -F 2048 \
  possorted_genome_bam.bam > possorted_genome_bam_filter.bam
samtools index possorted_genome_bam_filter.bam

# subsample 5% (NOTE: some lists mention ~2% / sub2.bam; command uses 5%)
samtools view -b -s 0.05 possorted_genome_bam_filter.bam > possorted_genome_bam_filter_sub5.bam
samtools index possorted_genome_bam_filter_sub5.bam

# --- 1.7 SMART-seq3 (split 5' UMI vs internal reads) ---
cd /path/to/floraseq/benchmark/smartseq3/zumi_fwdprimer_sub10_fig

samtools view -b -F 4 -F 256 -F 2048 \
  Smartseq3_fwdprimer.filtered.Aligned.GeneTagged.UBcorrected.sorted.bam \
  > Smartseq3_fwdprimer.filtered.Aligned.GeneTagged.UBcorrected.sorted.filter.bam
samtools index Smartseq3_fwdprimer.filtered.Aligned.GeneTagged.UBcorrected.sorted.filter.bam

# 5' UMI reads: UB tag matches exactly 8 letters
samtools view -h -e '[UB] =~ "^[A-Z]{8}$"' \
  Smartseq3_fwdprimer.filtered.Aligned.GeneTagged.UBcorrected.sorted.filter.bam \
  | samtools view -b - > Smartseq3.5prime_UMI.bam
samtools index Smartseq3.5prime_UMI.bam

# Internal reads: UB does not match 8-letter UMI
samtools view -h -e '[UB] !~ "^[A-Z]{8}$"' \
  Smartseq3_fwdprimer.filtered.Aligned.GeneTagged.UBcorrected.sorted.filter.bam \
  | samtools view -b - > Smartseq3.Internal.bam
samtools index Smartseq3.Internal.bam


# =============================================================================
# 2. RSeQC geneBody_coverage.py
#    Prepare one BAM path per line in cover_*.txt, then submit jobs.
# =============================================================================
#
# BAM list used for multi-tech comparison (cover_5tools.txt example):
#   .../mix3_human/outs/possorted_genome_bam_filter_sub10.bam
#   .../mixR_human/outs/possorted_genome_bam_filter.bam
#   .../vasaseq/star/vasa_293tAligned.sortedByCoord.out_sub5.bam
#   .../293t_mix_1k_v3-1_human/outs/possorted_genome_bam_filter_sub17.bam
#   .../10k_hgmm_5pv2_human/outs/possorted_genome_bam_filter_sub2.bam
#   .../zumi_fwdprimer_sub10_fig/Smartseq3.5prime_UMI.bam
#   .../zumi_fwdprimer_sub10_fig/Smartseq3.Internal.bam
#
# Reference BED:
#   /path/to/reference/RSeQC_ref/hg38_GENCODE.v38.bed

# --- 2.1 SMART-seq3 only (example SLURM job) ---
# cd .../zumi_fwdprimer_sub10_fig/coverage
# geneBody_coverage.py \
#   -r .../hg38_GENCODE.v38.bed \
#   -i cover_smartseq3.txt \
#   -o RseQC_smartseq3

# --- 2.2 Multi-tech comparison ---
# cd .../benchmark/coverage
# geneBody_coverage.py \
#   -r .../hg38_GENCODE.v38.bed \
#   -i cover_5tools.txt \
#   -o RseQC_5tools

# --- 2.3 10x-5 separately (if not included above) ---
# geneBody_coverage.py \
#   -r .../hg38_GENCODE.v38.bed \
#   -i cover_10x_5.txt \
#   -o RseQC_10x_5


# =============================================================================
# 3. Plot final benchmark figure
#    Kept: FLORA-seq-3, FLORA-seq-R, VASA-seq, Smartseq3.Internal
# =============================================================================
# python 2-plot_genebody.py
