#!/usr/bin/env bash
# =============================================================================
# Module : 3-mHypo / STARsolo Gene + SJ (group1..24_3r)
# Source : 3_mHypo/bin/old/MARVEL_DropBased.sh (STARsolo section)
# Next: Rscript 8-sj-marvel.R  (needs meta.csv from 4-anno-coarse.py)
# Environment: STAR, samtools
# =============================================================================

set -euo pipefail

WHITELIST="/path/to/cellranger/lib/python/cellranger/barcodes/737K-august-2016.txt"
GENOME_DIR="/path/to/floraseq/marvel/basic/gtf/starsolo_mm_index"
WORK_DIR="/path/to/floraseq/marvel/proj_mHypo"
BAM_ROOT="${WORK_DIR}/cellranger"

cd "${WORK_DIR}"
mkdir -p ./star

# Original loop ran groups 9..24; full 1..24 recommended for a complete re-run.
for i in $(seq 1 24); do
  sample="group${i}_3r"
  echo "Running STARsolo: ${sample}"
  STAR --runThreadN 16 \
       --genomeDir "${GENOME_DIR}" \
       --soloType CB_UMI_Simple \
       --readFilesIn "${BAM_ROOT}/${sample}/outs/possorted_genome_bam.bam" \
       --readFilesCommand "samtools view -F 0x100" \
       --readFilesType SAM SE \
       --soloInputSAMattrBarcodeSeq CR UR \
       --soloInputSAMattrBarcodeQual CY UY \
       --soloCBwhitelist "${WHITELIST}" \
       --soloFeatures Gene SJ

  mkdir -p "./star/${sample}"
  rm -f Aligned.out.sam
  mv Log.final.out Log.out Log.progress.out SJ.out.tab ./Solo.out
  mv Solo.out "./star/${sample}"
done

echo "STARsolo done. Next: Rscript 8-sj-marvel.R"
