#!/usr/bin/env bash
# Copy this file to config.sh, edit the paths, then run: source config.sh

export FLORASEQ_DATA_ROOT="/path/to/floraseq-data"
export FLORASEQ_REFERENCE_ROOT="/path/to/reference-files"
export FLORASEQ_RESULTS_ROOT="/path/to/floraseq-results"

# Optional overrides used by the five-method sensitivity workflow.
# export Q6_SMART_BAM="/path/to/smartseq3.bam"
# export Q6_SMART_WHITELIST="/path/to/smartseq3_barcodes.txt"
# export Q6_VASA_R1="/path/to/VASA_R1.fastq.gz"
# export Q6_VASA_R2="/path/to/VASA_R2.fastq.gz"
