#!/usr/bin/env bash
# =============================================================================
# Module : 0-download / Public benchmark data download
# Dataset: Cell-line benchmark (VASA-seq, 10x 3'/5', FLASH-seq notes)
# Source : 1_Cell_line/0-bin/1-down-data.sh
# =============================================================================

# =============================================================================
# 1. VASA-seq
# =============================================================================
# Step 1: download .sra
cd /path/to/floraseq/benchmark/vasaseq/fq/
prefetch SRR14783059

# Step 2: convert .sra to .fastq
fasterq-dump --split-files SRR14783059 -e 16

# Step 3: compress (recommended)
gzip SRR14783059*.fastq


# =============================================================================
# 2. 10x 3'
# =============================================================================
# 293T + NIH3T3 (1k cell) - 3' v3.1:
wget https://www.10xgenomics.com/datasets/1-k-1-1-mixture-of-human-hek-293-t-and-mouse-nih-3-t-3-cells-3-v-3-1-3-1-standard-6-0-0


# =============================================================================
# 3. 10x 5'
# =============================================================================
# 293T + NIH3T3:
# https://www.10xgenomics.com/datasets/10-k-1-1-mixture-of-human-hek-293-t-and-mouse-nih-3-t-3-cells-5-v-2-0-chromium-x-2-standard-6-1-0
# NOTE [KNOWN_ISSUES]: the wget target below points to a dataset landing page,
# not a direct FASTQ archive URL. Replace with the real download link if needed.
wget https://www.10xgenomics.com/datasets/10-k-1-1-mixture-of-human-hek-293-t-and-mouse-nih-3-t-3-cells-5-v-2-0-chromium-x-2-standard-6-1-0

# =============================================================================
# 4. SMART-seq3
# =============================================================================

