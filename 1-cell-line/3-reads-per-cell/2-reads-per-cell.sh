#!/usr/bin/env bash
# =============================================================================
# Module : 3-reads-per-cell / Sequencing saturation (HEK293T, human BAM)
# Dataset: FLORA-seq, 10x-3', 10x-5' (final benchmark figure)
# Source : 1_Cell_line/0-bin/4-reads-per-cell-human-only.sh
#
# Workflow:
#   1) subset-bam → HEK293T-only BAM (barcodes from 1-get-barcode.py / 5-matrix)
#   2) Downsample along a mean-reads-per-cell gradient
#   3) Count unique UMI / genes per cell (CB / UB / GX)
#   4) Plot: python 3-plot_reads_per_cell.py
#
# Usage (run saturation once per method):
#   bash 2-reads-per-cell.sh FLORA-seq
#   bash 2-reads-per-cell.sh 10x-3
#   bash 2-reads-per-cell.sh 10x-5
# =============================================================================

set -euo pipefail

METHOD="${1:-}"

# Depth gradient (mean reads per cell)
DEPTHS=(2500 5000 10000 20000 30000 40000 50000 60000)

# =============================================================================
# 1. subset-bam: HEK293T-only BAMs
#    Run once:  bash 2-reads-per-cell.sh
# =============================================================================
# Barcode lists (from 1-get-barcode.py / 5-matrix/3-process-matrix.py):
#   FLORA: .../mixR3/outs/hek293t_barcode.txt
#   10x-3: .../293t_mix_1k_v3-1_run/outs/hek293t_barcode.txt
#   10x-5: .../10k_hgmm_5pv2_run/outs/hek293t_barcode.txt

if [[ -z "${METHOD}" ]]; then
  # FLORA-seq (human-mapped BAM)
  cd /path/to/floraseq/cell_line/mixR3_human/outs/
  subset-bam --bam possorted_genome_bam.bam \
    --cell-barcodes /path/to/floraseq/cell_line/mixR3/outs/hek293t_barcode.txt \
    --out-bam flora_hek293t.bam

  # 10x-3'
  cd /path/to/floraseq/benchmark/10x_3/293t_mix_1k_v3-1/293t_mix_1k_v3-1_human/outs
  subset-bam --bam possorted_genome_bam.bam \
    --cell-barcodes /path/to/floraseq/benchmark/10x_3/293t_mix_1k_v3-1/293t_mix_1k_v3-1_run/outs/hek293t_barcode.txt \
    --out-bam 10x_3_hek293t.bam

  # 10x-5'
  cd /path/to/floraseq/benchmark/10x_5/10k_hgmm_5pv2_human/outs
  subset-bam --bam possorted_genome_bam.bam \
    --cell-barcodes /path/to/floraseq/benchmark/10x_5/10k_hgmm_5pv2_run/outs/hek293t_barcode.txt \
    --out-bam 10x_5_hek293t.bam

  echo "subset-bam done."
  echo "Next: bash 2-reads-per-cell.sh {FLORA-seq|10x-3|10x-5}"
  exit 0
fi


# =============================================================================
# 2. Saturation loop (one method per run)
# =============================================================================
# N_CELLS = barcode count used to convert target mean depth → total reads.

case "${METHOD}" in
  FLORA-seq)
    mkdir -p /path/to/floraseq/cell_line/mixR3_human/outs/depth
    cd /path/to/floraseq/cell_line/mixR3_human/outs/depth
    INPUT_BAM="/path/to/floraseq/cell_line/mixR3_human/outs/flora_hek293t.bam"
    METHOD_NAME="FLORA-seq"
    OUTPUT_CSV="/path/to/floraseq/cell_line/mixR3_human/outs/depth/flora_stats.csv"
    N_CELLS=351
    ;;
  10x-3)
    mkdir -p /path/to/floraseq/benchmark/10x_3/293t_mix_1k_v3-1/293t_mix_1k_v3-1_human/outs/depth
    cd /path/to/floraseq/benchmark/10x_3/293t_mix_1k_v3-1/293t_mix_1k_v3-1_human/outs/depth
    INPUT_BAM="/path/to/floraseq/benchmark/10x_3/293t_mix_1k_v3-1/293t_mix_1k_v3-1_human/outs/10x_3_hek293t.bam"
    METHOD_NAME="10x-3"
    OUTPUT_CSV="/path/to/floraseq/benchmark/10x_3/293t_mix_1k_v3-1/293t_mix_1k_v3-1_human/outs/depth/10x_3_stats.csv"
    N_CELLS=505
    ;;
  10x-5)
    mkdir -p /path/to/floraseq/benchmark/10x_5/10k_hgmm_5pv2_human/outs/depth
    cd /path/to/floraseq/benchmark/10x_5/10k_hgmm_5pv2_human/outs/depth
    INPUT_BAM="/path/to/floraseq/benchmark/10x_5/10k_hgmm_5pv2_human/outs/10x_5_hek293t.bam"
    METHOD_NAME="10x-5"
    OUTPUT_CSV="/path/to/floraseq/benchmark/10x_5/10k_hgmm_5pv2_human/outs/depth/10x_5_stats.csv"
    N_CELLS=5564
    ;;
  *)
    echo "Unknown METHOD: ${METHOD}"
    echo "Usage:"
    echo "  bash 2-reads-per-cell.sh                 # subset-bam only"
    echo "  bash 2-reads-per-cell.sh FLORA-seq|10x-3|10x-5"
    exit 1
    ;;
esac

echo "Method,TargetDepth,CellID,UMI_Count,Gene_Count" > "${OUTPUT_CSV}"

echo "Calculating total reads for ${METHOD_NAME}..."
TOTAL_READS=$(samtools view -c "${INPUT_BAM}")
echo "Total Reads: ${TOTAL_READS} | N_CELLS=${N_CELLS}"

for DEPTH in "${DEPTHS[@]}"; do
    TARGET_TOTAL=$((DEPTH * N_CELLS))
    FACTOR=$(awk -v target="${TARGET_TOTAL}" -v total="${TOTAL_READS}" 'BEGIN {print target / total}')

    if (( $(echo "${FACTOR} > 1.0" | bc -l) )); then
        echo "Skipping depth ${DEPTH}: not enough reads (need ${TARGET_TOTAL}, have ${TOTAL_READS})"
        continue
    fi

    # samtools -s: integer seed + fraction (seed=42)
    SUBSAMPLE_PARAM="42${FACTOR}"
    echo "Processing depth=${DEPTH} (factor=${FACTOR})..."

    # Stream: subsample → keep CB/UB/GX → unique UMI per cell-gene → counts
    samtools view -s "${SUBSAMPLE_PARAM}" "${INPUT_BAM}" \
    | grep "CB:Z:" | grep "UB:Z:" | grep "GX:Z:" \
    | awk '{
        cb=""; ub=""; gx="";
        for (i = 12; i <= NF; i++) {
            if ($i ~ /^CB:Z:/) cb = substr($i, 6);
            if ($i ~ /^UB:Z:/) ub = substr($i, 6);
            if ($i ~ /^GX:Z:/) gx = substr($i, 6);
        }
        if (cb != "" && ub != "" && gx != "") print cb, gx, ub
    }' \
    | sort -k1,1 -k2,2 -k3,3 -u \
    | awk -v method="${METHOD_NAME}" -v depth="${DEPTH}" '
    {
        c = $1; g = $2;
        umi_count[c]++;
        if (!((c SUBSEP g) in seen_gene)) {
            gene_count[c]++;
            seen_gene[c, g] = 1;
        }
    }
    END {
        for (c in umi_count) {
            print method "," depth "," c "," umi_count[c] "," gene_count[c]
        }
    }' >> "${OUTPUT_CSV}"
done

echo "Done! Results → ${OUTPUT_CSV}"
echo "After all methods finish: python 3-plot_reads_per_cell.py"
