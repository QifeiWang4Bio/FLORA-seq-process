#!/usr/bin/env bash
# =============================================================================
# Module : 4-sj / FLORA-seq vs long-read splice-junction preparation
# Dataset: FLORA-seq K562 (mixR3 / mix3 / mixR) + ENCODE / SGNex PacBio BAMs
# Source : shell blocks in 1_Cell_line/0-bin/14-long-read-benchmark.py
#
# Workflow:
#   1) Stage long-read BAMs (ENCODE K562; see notes below)
#   2) subset-bam FLORA-seq BAMs to K562 barcodes
#   3) Extract splice junctions with regtools
#   4) Validate with 4-long-read-sj-validation.py
#
# Final figure (from validation script):
#   .../benchmark/long_read/fig/Fig_LongRead_Validation_Panel.pdf
# =============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BARCODES="/path/to/floraseq/cell_line/mixR3/outs/k562_barcode.csv"
LONG_DIR="/path/to/floraseq/benchmark/long_read"
FIG_DIR="${LONG_DIR}/fig"

REGTOOLS_OPTS=(-s XS -a 8 -m 50 -M 500000)
# -s XS : strand from XS tag
# -a 8  : minimum anchor length (bp)
# -m/-M : intron length range


# =============================================================================
# 1. Long-read data (ENCODE K562) — stage under LONG_DIR
# =============================================================================
# ENCFF322UJU  ENCSR526TQU  k562_1_1
# ENCFF504GVG  ENCSR589FUJ  k562_2_1
# ENCFF479SQR  ENCSR589FUJ  k562_2_2
# ENCFF661OEY  ENCSR983KDL  k562_3_1
# ENCFF645UVN  ENCSR983KDL  k562_3_2
# https://www.encodeproject.org/experiments/ENCSR526TQU/#browser


# =============================================================================
# 2. Subset FLORA-seq BAMs to K562 barcodes
# =============================================================================

# 3'+Random
cd /path/to/floraseq/cell_line/mixR3_human/outs/
subset-bam --bam possorted_genome_bam.bam \
  --cell-barcodes "${BARCODES}" \
  --out-bam flora_k562.bam
samtools index flora_k562.bam

# Random only
cd /path/to/floraseq/cell_line/mixR_human/outs/
subset-bam --bam possorted_genome_bam.bam \
  --cell-barcodes "${BARCODES}" \
  --out-bam flora_k562.bam
samtools index flora_k562.bam

# 3' only
cd /path/to/floraseq/cell_line/mix3_human/outs/
subset-bam --bam possorted_genome_bam.bam \
  --cell-barcodes "${BARCODES}" \
  --out-bam flora_k562.bam
samtools index flora_k562.bam


# =============================================================================
# 3. Extract splice junctions (regtools)
# =============================================================================

# --- 3.1 FLORA-seq ---
cd /path/to/floraseq/cell_line/mixR3_human/outs/
regtools junctions extract "${REGTOOLS_OPTS[@]}" \
  flora_k562.bam > flora_k562_sjs.bed
# expected ~136946 junctions for 3'+Random

cd /path/to/floraseq/cell_line/mix3_human/outs/
regtools junctions extract "${REGTOOLS_OPTS[@]}" \
  flora_k562.bam > flora_k562_sjs.bed

cd /path/to/floraseq/cell_line/mixR_human/outs/
regtools junctions extract "${REGTOOLS_OPTS[@]}" \
  flora_k562.bam > flora_k562_sjs.bed

# --- 3.2 ENCODE long-read BAMs ---
cd "${LONG_DIR}"

bam_files=(
  "ENCFF322UJU.bam"
  "ENCFF479SQR.bam"
  "ENCFF504GVG.bam"
  "ENCFF645UVN.bam"
  "ENCFF661OEY.bam"
)

for bam_file in "${bam_files[@]}"; do
  bed_file="${bam_file%.bam}_sjs.bed"
  echo "Processing: ${bam_file} -> ${bed_file}"
  samtools index "${bam_file}"
  regtools junctions extract "${REGTOOLS_OPTS[@]}" \
    "${bam_file}" > "${bed_file}"
done

# --- 3.3 SGNex PacBio (optional alternative reference; not used in final plot) ---
# bam_files=(
#   "SGNex_K562_PacBio-SMRTcell_replicate7_run1.bam"
# )
# for bam_file in "${bam_files[@]}"; do
#   bed_file="${bam_file%.bam}_sjs.bed"
#   echo "Processing: ${bam_file} -> ${bed_file}"
#   samtools index "${bam_file}"
#   regtools junctions extract "${REGTOOLS_OPTS[@]}" \
#     "${bam_file}" > "${bed_file}"
# done


# =============================================================================
# 4. Validate FLORA SJs against long-read SJs
# =============================================================================
# Default FLORA_BED in the Python script is mixR3 (3'+Random).
# Switch FLORA_BED there to compare Random-only or 3'-only libraries.
mkdir -p "${FIG_DIR}"
python3 "${SCRIPT_DIR}/4-long-read-sj-validation.py"

echo "Done. Final figure: ${FIG_DIR}/Fig_LongRead_Validation_Panel.pdf"
