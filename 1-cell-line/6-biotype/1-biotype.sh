#!/usr/bin/env bash
# =============================================================================
# Module : 6-biotype / Cross-technology gene biotype fractions (RCAS)
# Dataset: Cell-line benchmark (human-mapped, depth-matched BAMs)
# Source : 1_Cell_line/0-bin/8-biotype.sh (+ 8-biotype-rcas.R)
#
# Workflow:
#   1) BAM → BED (bedtools bamtobed)
#   2) RCAS biotype summary (2-biotype-rcas.R)
#
# BAMs match the genebody subsample set (comparable read depth).
# GTF: Ensembl GRCh38.109 with Mt gene annotations (gene_biotype).
#
# Outputs (per sample, under <bam_dir>/biotype/):
#   RCAS_<name>.csv
#   RCAS_<name>.pdf
#
# Usage:
#   bash 1-biotype.sh
#   # or call R directly after BED files exist:
#   Rscript 2-biotype-rcas.R "$GTF_FILE" bed1.bed bed2.bed ...
# =============================================================================

set -euo pipefail

# Directory of this script (for 2-biotype-rcas.R)
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

GTF_FILE="/path/to/reference/ensembl/Homo_sapiens.GRCh38.109.chr_Mtgene.gtf"

# Depth-matched BAMs (same set as 2-genebody)
BAM_FILES=(
  "/path/to/floraseq/cell_line/mixR3_human/outs/possorted_genome_bam_filter_sub7.bam"
  "/path/to/floraseq/cell_line/mix3_human/outs/possorted_genome_bam_filter_sub10.bam"
  "/path/to/floraseq/cell_line/mixR_human/outs/possorted_genome_bam_filter.bam"
  "/path/to/floraseq/benchmark/vasaseq/star/vasa_293tAligned.sortedByCoord.out_sub5.bam"
  "/path/to/floraseq/benchmark/10x_3/293t_mix_1k_v3-1/293t_mix_1k_v3-1_human/outs/possorted_genome_bam_filter_sub17.bam"
  "/path/to/floraseq/benchmark/10x_5/10k_hgmm_5pv2_human/outs/possorted_genome_bam_filter_sub2.bam"
  "/path/to/floraseq/benchmark/smartseq3/zumi_fwdprimer_sub10_fig/Smartseq3_fwdprimer.filtered.Aligned.GeneTagged.UBcorrected.sorted.filter.bam"
)


# =============================================================================
# 1. BAM → BED
# =============================================================================
BED_FILES=()
for bam_file in "${BAM_FILES[@]}"; do
  bed_file="${bam_file%.bam}.bed"
  BED_FILES+=("${bed_file}")

  if [[ -f "${bed_file}" ]]; then
    echo "BED exists, skip: ${bed_file}"
    continue
  fi

  echo "bamtobed: ${bam_file}"
  bedtools bamtobed -i "${bam_file}" > "${bed_file}"
  echo "Wrote: ${bed_file}"
done

echo "BAM → BED done (${#BED_FILES[@]} files)."


# =============================================================================
# 2. RCAS biotype fraction
# =============================================================================
echo "Running RCAS with GTF: ${GTF_FILE}"
Rscript "${SCRIPT_DIR}/2-biotype-rcas.R" "${GTF_FILE}" "${BED_FILES[@]}"

echo "All biotype analyses completed."
