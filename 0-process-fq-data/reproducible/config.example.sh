#!/usr/bin/env bash

# Copy this file, replace the placeholder paths, and pass it to a runner.
export FLORA_FASTQ_OUTPUT_ROOT="/path/to/output"

export CELL_LINE_RAW_R1="/path/to/scMix_3R_L01_R1.fq.gz"
export CELL_LINE_RAW_R2="/path/to/scMix_3R_L01_R2.fq.gz"

export CAF_BATCH1_RAW_R1="/path/to/2P221113028US2S2718DX_scB_P3_L01_R1.fq.gz"
export CAF_BATCH1_RAW_R2="/path/to/2P221113028US2S2718DX_scB_P3_L01_R2.fq.gz"
export CAF_BATCH2_RAW_R1="/path/to/2P221114048US2S2699DX_scB_P3_L01_R1.fq.gz"
export CAF_BATCH2_RAW_R2="/path/to/2P221114048US2S2699DX_scB_P3_L01_R2.fq.gz"

export MHYPO_BATCH1_RAW_R1="/path/to/snHy1_R1.fq.gz"
export MHYPO_BATCH1_RAW_R2="/path/to/snHy1_R2.fq.gz"
export MHYPO_BATCH2_RAW_R1="/path/to/snHy2_R1.fq.gz"
export MHYPO_BATCH2_RAW_R2="/path/to/snHy2_R2.fq.gz"

# Executable names can be replaced by full paths if they are not on PATH.
export FASTP_BIN="fastp"

# Cell-line options recovered from the archived run script.
CELL_LINE_FASTP_ARGS=(-q 20 -u 30 -n 10 -w 4)
CELL_LINE_RANDOM_FASTP_ARGS=(-q 20 -u 30 -n 10 -w 4 -l 150)
# "keep" matches production Cleandata (mixR R1 = 158/142 bp, i.e. only the RT
# index is removed, the downstream sequence is not truncated). Set to "trim"
# only to reproduce a hypothetical fixed-length variant for comparison.
CELL_LINE_RANDOM_TAIL_MODE="keep"
# These settings are explicitly recorded in both mHypo QC scripts.
MHYPO_FASTP_ARGS=(-q 20 -u 30 -n 10 -w 12 -l 26 --cut_right)
export RUN_FASTQC=1
