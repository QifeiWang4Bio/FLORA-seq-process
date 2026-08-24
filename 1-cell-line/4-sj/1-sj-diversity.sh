#!/usr/bin/env bash
# =============================================================================
# Module : 4-sj / Junction saturation (SJ diversity) across technologies
# Source : shell blocks in 1_Cell_line/0-bin/13-transcipt-diverisity.py
#
# Workflow:
#   1) Depth-match BAMs to ~7.3–7.5 M usable reads (where needed)
#   2) Run RSeQC junction_saturation.py
#   3) Plot with 2-plot_sj_diversity.py
#
# Final figures (working dir .../benchmark/sj_diversity/):
#   Benchmark_Junctions_All.pdf
#   Benchmark_Junctions_Known.pdf
#   Benchmark_Junctions_Novel.pdf
#
# Methods in final plot: FLORA-seq (3'), FLORA-seq (Random), 10x-3', Smart-seq3
# (VASA saturation is also computed below but not included in the final plot.)
# =============================================================================

set -euo pipefail

REF_BED="/path/to/reference/RSeQC_ref/hg38_GENCODE.v38.bed"
OUT_DIR="/path/to/floraseq/benchmark/sj_diversity"
mkdir -p "${OUT_DIR}"


# =============================================================================
# 1. Depth-match BAMs (target ~7.3–7.5 M usable reads)
# =============================================================================
# Approximate usable-read counts before this step:
#   FLORA-seq Random : .../mixR_human/.../possorted_genome_bam_filter.bam           ~7,421,528
#   FLORA-seq 3'     : .../mix3_human/.../possorted_genome_bam_filter_sub10.bam     ~10,210,256
#   10x-3'           : .../possorted_genome_bam_filter_sub17.bam                    ~7,504,513
#   VASA-seq         : .../vasa_293tAligned.sortedByCoord.out_sub5.bam              ~7,262,441
#   SMART-seq3       : .../zumi_fwdprimer_sub10_fig/...sorted.filter.bam            ~21,628,873
#
# Only FLORA-seq 3' and SMART-seq3 need an extra samtools subsample here.
# Seed prefix "42." is the samtools -s random seed.

# --- FLORA-seq 3' : 10,210,256 → ~7.4 M ---
cd /path/to/floraseq/cell_line/mix3_human/outs/
samtools view -@ 8 -b -s 42.727 \
  possorted_genome_bam_filter_sub10.bam \
  > possorted_genome_bam_for_sjDiversity.bam
samtools index possorted_genome_bam_for_sjDiversity.bam

# --- SMART-seq3 : 21,628,873 → ~7.4 M ---
cd /path/to/floraseq/benchmark/smartseq3/zumi_fwdprimer_sub10_fig/
samtools view -@ 8 -b -s 42.343 \
  Smartseq3_fwdprimer.filtered.Aligned.GeneTagged.UBcorrected.sorted.filter.bam \
  > possorted_genome_bam_for_sjDiversity.bam
samtools index possorted_genome_bam_for_sjDiversity.bam


# =============================================================================
# 2. junction_saturation.py (RSeQC)
# =============================================================================
cd "${OUT_DIR}"

# 1) FLORA-seq Random (already ~7.4 M; no extra subsample)
junction_saturation.py \
  -i /path/to/floraseq/cell_line/mixR_human/outs/possorted_genome_bam_filter.bam \
  -r "${REF_BED}" \
  -o flora_human_R_saturation

# 2) FLORA-seq 3'
junction_saturation.py \
  -i /path/to/floraseq/cell_line/mix3_human/outs/possorted_genome_bam_for_sjDiversity.bam \
  -r "${REF_BED}" \
  -o flora_human_3_saturation

# 3) 10x-3'
junction_saturation.py \
  -i /path/to/floraseq/benchmark/10x_3/293t_mix_1k_v3-1/293t_mix_1k_v3-1_human/outs/possorted_genome_bam_filter_sub17.bam \
  -r "${REF_BED}" \
  -o 10x_3_human_saturation

# 4) VASA-seq (saturation computed; not in final plot)
junction_saturation.py \
  -i /path/to/floraseq/benchmark/vasaseq/star/vasa_293tAligned.sortedByCoord.out_sub5.bam \
  -r "${REF_BED}" \
  -o vasa_293t_saturation

# 5) SMART-seq3
junction_saturation.py \
  -i /path/to/floraseq/benchmark/smartseq3/zumi_fwdprimer_sub10_fig/possorted_genome_bam_for_sjDiversity.bam \
  -r "${REF_BED}" \
  -o smartseq3_saturation


# =============================================================================
# 3. Plot
# =============================================================================
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
python3 "${SCRIPT_DIR}/2-plot_sj_diversity.py"

echo "Done. Final figures under ${OUT_DIR}/:"
echo "  Benchmark_Junctions_All.pdf"
echo "  Benchmark_Junctions_Known.pdf"
echo "  Benchmark_Junctions_Novel.pdf"
