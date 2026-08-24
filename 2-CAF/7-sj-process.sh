#!/usr/bin/env bash
# =============================================================================
# Module : 2-CAF / STARsolo Gene + SJ quantification for CAF SJ analysis
# Dataset: Breast cancer FLORA-seq (human_sample1/2/3; library designs 3r / 3 / r)
# Source : 2_hBreast/bin/5-sj-process.sh
#
# This is step 4 prep of the CAF track (before MARVEL):
#   Input  : Cell Ranger possorted_genome_bam.bam per sample
#   Output : STARsolo Solo.out/Gene + Solo.out/SJ under caf_sc_update/
#
# Final path used by 8-sj-marvel.R:
#   human_sample{1,2,3}_3r   (3'+Random merged library)
#
# Optional tracks (_3 polyA-only, _r random-only) kept below for completeness.
#
# Prerequisite:
#   CAF barcodes / subtypes → meta.csv from 6-caf-analysis.py
#   (MARVEL subsets cells with that metadata; STARsolo itself uses CB whitelist)
# =============================================================================

set -euo pipefail

source /path/to/software/miniforge3/etc/profile.d/conda.sh
conda activate sc-bio

WHITELIST="/path/to/software/cellranger/lib/python/cellranger/barcodes/737K-august-2016.txt"
GENOME_DIR="/path/to/reference/starsolo_human_index"
WORK_DIR="/path/to/floraseq/marvel/caf_sc_update"
BAM_ROOT="/path/to/floraseq/caf"

# Reference notes (from original notebook):
#   GTF: .../marvel/basic/gtf/gencode.v31.annotation.gtf
#   Cell Ranger whitelist: 737K-august-2016.txt
#   STAR genomeGenerate (if rebuilding index):
#     STAR --runMode genomeGenerate --runThreadN 24 \
#       --genomeDir ./starsolo_human_index \
#       --genomeFastaFiles GRCh38.primary_assembly.genome.fa \
#       --sjdbGTFfile gencode.v31.annotation.gtf

run_starsolo() {
  local sample=$1
  local threads=${2:-20}

  echo "Running STARsolo: ${sample} ($(date))"
  mkdir -p "${WORK_DIR}/${sample}"
  cd "${WORK_DIR}/${sample}"

  STAR --runThreadN "${threads}" \
    --genomeDir "${GENOME_DIR}" \
    --soloType CB_UMI_Simple \
    --readFilesIn "${BAM_ROOT}/${sample}/outs/possorted_genome_bam.bam" \
    --readFilesCommand samtools view -F 0x100 \
    --readFilesType SAM SE \
    --soloInputSAMattrBarcodeSeq CR UR \
    --soloInputSAMattrBarcodeQual CY UY \
    --soloCBwhitelist "${WHITELIST}" \
    --soloFeatures Gene SJ

  rm -f Aligned.out.sam
}


# =============================================================================
# 1. FINAL — 3'+Random (human_sample*_3r)  → used by 8-sj-marvel.R
# =============================================================================
# SLURM example (optional):
#   #SBATCH --cpus-per-task=20
#   #SBATCH --mem=120G
#   #SBATCH --partition=corexd192

cd "${WORK_DIR}"
for i in 1 2 3; do
  run_starsolo "human_sample${i}_3r" 20
done
echo "3r STARsolo finished on $(date)"


# =============================================================================
# 2. OPTIONAL — polyA-only (human_sample*_3)
# =============================================================================
# for i in 1 2 3; do
#   run_starsolo "human_sample${i}_3" 20
# done


# =============================================================================
# 3. OPTIONAL — random-only (human_sample*_r)
# =============================================================================
# for i in 1 2 3; do
#   run_starsolo "human_sample${i}_r" 12
# done

echo "Done. Next: Rscript 8-sj-marvel.R  (expects Solo.out under ${WORK_DIR})"
