#!/usr/bin/env bash
# =============================================================================
# Module : 6-biotype / Per-cell biotype fractions (FLORA-seq vs 10x-3)
# Source : 1_Cell_line/0-bin/8-biotype-per-cell.sh
#
# Prerequisite:
#   HEK293T-only BAMs from 3-reads-per-cell/2-reads-per-cell.sh
#     flora_hek293t.bam
#     10x_3_hek293t.bam
#
# Workflow:
#   1) Index BAMs
#   2) BAM → BED with CB barcode in column 4
#   3) RCAS annotate once per dataset (4-biotype-percell-rcas.R)
#   4) Plot per-cell histograms (5-plot_biotype_percell.py)
#
# Final figure:
#   .../mixR3_human/outs/biotype_percell/biotype_percell_hist_median.pdf
# =============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

GTF="/path/to/reference/ensembl/Homo_sapiens.GRCh38.109.chr_Mtgene.gtf"

FLORA_BAM="/path/to/floraseq/cell_line/mixR3_human/outs/flora_hek293t.bam"
FLORA_OUTDIR="/path/to/floraseq/cell_line/mixR3_human/outs/biotype_percell"

TENX_BAM="/path/to/floraseq/benchmark/10x_3/293t_mix_1k_v3-1/293t_mix_1k_v3-1_human/outs/10x_3_hek293t.bam"
TENX_OUTDIR="/path/to/floraseq/benchmark/10x_3/293t_mix_1k_v3-1/293t_mix_1k_v3-1_human/outs/biotype_percell"

mkdir -p "${FLORA_OUTDIR}" "${TENX_OUTDIR}"


# =============================================================================
# 1. Index HEK293T BAMs
# =============================================================================
samtools index "${FLORA_BAM}"
samtools index "${TENX_BAM}"


# =============================================================================
# 2. BAM → BED  (CB in col 4; CIGAR-derived reference span; strand from FLAG)
# =============================================================================
# Mapped reads only (-F 4); skip reads without CB:Z:.
# ref_span sums CIGAR ops that consume reference (M D N = X).

bam_to_cb_bed() {
    local BAM=$1 OUT=$2
    samtools view -F 4 "${BAM}" | gawk '
    function ref_span(cig,    len) {
        len = 0
        while (match(cig, /[0-9]+[A-Z]/)) {
            op = substr(cig, RSTART + RLENGTH - 1, 1)
            if (op ~ /[MDN=X]/) len += substr(cig, RSTART, RLENGTH - 1)
            cig = substr(cig, RSTART + RLENGTH)
        }
        return len
    }
    {
        cb = ""
        for (i = 12; i <= NF; i++)
            if (substr($i, 1, 5) == "CB:Z:") { cb = substr($i, 6); break }
        if (cb == "") next
        start  = $4 - 1
        end    = start + ref_span($6)
        strand = (and($2, 16)) ? "-" : "+"
        print $3 "\t" start "\t" end "\t" cb "\t0\t" strand
    }' > "${OUT}"
    printf "[%s] BED written: %s  (%d reads)\n" \
        "$(date +%H:%M:%S)" "${OUT}" "$(wc -l < "${OUT}")"
}

bam_to_cb_bed "${FLORA_BAM}" "${FLORA_OUTDIR}/flora_cb.bed"
bam_to_cb_bed "${TENX_BAM}"  "${TENX_OUTDIR}/10x3_cb.bed"


# =============================================================================
# 3. Per-read biotype annotation (RCAS, once per dataset)
# =============================================================================
Rscript "${SCRIPT_DIR}/4-biotype-percell-rcas.R" \
    "${GTF}" "${FLORA_OUTDIR}/flora_cb.bed" "${FLORA_OUTDIR}/flora_read_biotype.csv"

Rscript "${SCRIPT_DIR}/4-biotype-percell-rcas.R" \
    "${GTF}" "${TENX_OUTDIR}/10x3_cb.bed" "${TENX_OUTDIR}/10x3_read_biotype.csv"


# =============================================================================
# 4. Per-cell fraction histograms
#    Final figure: biotype_percell_hist_median.pdf
# =============================================================================
python3 "${SCRIPT_DIR}/5-plot_biotype_percell.py"

echo "Done. Final figure:"
echo "  ${FLORA_OUTDIR}/biotype_percell_hist_median.pdf"
