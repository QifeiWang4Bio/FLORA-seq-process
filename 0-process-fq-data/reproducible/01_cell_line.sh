#!/usr/bin/env bash

set -euo pipefail
SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
source "${SCRIPT_DIR}/common.sh"

if [[ $# -ne 1 ]]; then
    echo "Usage: bash $0 /path/to/config.sh" >&2
    exit 1
fi
source "$1"

: "${FLORA_FASTQ_OUTPUT_ROOT:?Set FLORA_FASTQ_OUTPUT_ROOT in the config file}"
: "${CELL_LINE_RAW_R1:?Set CELL_LINE_RAW_R1 in the config file}"
: "${CELL_LINE_RAW_R2:?Set CELL_LINE_RAW_R2 in the config file}"
if ! declare -p CELL_LINE_FASTP_ARGS >/dev/null 2>&1; then
    CELL_LINE_FASTP_ARGS=(-q 20 -u 30 -n 10 -w 4)
fi
if ! declare -p CELL_LINE_RANDOM_FASTP_ARGS >/dev/null 2>&1; then
    CELL_LINE_RANDOM_FASTP_ARGS=(-q 20 -u 30 -n 10 -w 4 -l 150)
fi
# The archived pipeline only strips the RT index from the random-primer
# barcode read (Transfor_barcode.py + cut_index_R.py); it never truncates the
# downstream sequence away. Verified against production Cleandata: mixR R1 is
# 158/142 bp, not 26 bp. Keep this configurable in case a fixed-length
# variant is ever needed for comparison.
CELL_LINE_RANDOM_TAIL_MODE=${CELL_LINE_RANDOM_TAIL_MODE:-keep}
fastp_command=${FASTP_BIN:-fastp}

require_file "${CELL_LINE_RAW_R1}"
require_file "${CELL_LINE_RAW_R2}"
require_command python3
require_command "${fastp_command}"

work_dir="${FLORA_FASTQ_OUTPUT_ROOT}/cell_line/work"
final_dir="${FLORA_FASTQ_OUTPUT_ROOT}/cell_line/cellranger_fastq"
mkdir -p "${work_dir}" "${final_dir}"

python3 "${REFORMAT_PY}" reformat \
    --raw-r1 "${CELL_LINE_RAW_R1}" \
    --raw-r2 "${CELL_LINE_RAW_R2}" \
    --indices "${SCRIPT_DIR}/configs/cell_line_indices.tsv" \
    --outdir "${work_dir}/reformatted" \
    --tail-mode trim --gzip-output

three_r1="${work_dir}/reformatted/mix_3/mix_3_S1_L001_R1_001.fastq.gz"
three_r2="${work_dir}/reformatted/mix_3/mix_3_S1_L001_R2_001.fastq.gz"
random_r1="${work_dir}/reformatted/mix_r/mix_r_S1_L001_R1_001.fastq.gz"
random_r2="${work_dir}/reformatted/mix_r/mix_r_S1_L001_R2_001.fastq.gz"

python3 "${REFORMAT_PY}" split-raw \
    --raw-r1 "${CELL_LINE_RAW_R1}" \
    --raw-r2 "${CELL_LINE_RAW_R2}" \
    --indices "${SCRIPT_DIR}/configs/cell_line_indices.tsv" \
    --only-primer 3 \
    --outdir "${work_dir}/three_prime_raw" --gzip-output

three_raw_r1="${work_dir}/three_prime_raw/mix_3/mix_3_sequencer_R1.fastq.gz"
three_raw_r2="${work_dir}/three_prime_raw/mix_3/mix_3_sequencer_R2.fastq.gz"
python3 "${REFORMAT_PY}" reformat \
    --raw-r1 "${three_raw_r1}" \
    --raw-r2 "${three_raw_r2}" \
    --indices "${SCRIPT_DIR}/configs/cell_line_indices.tsv" \
    --only-primer 3 --cell-barcode-orientation raw \
    --outdir "${work_dir}/three_prime_pre_fastp" \
    --tail-mode trim --gzip-output

three_pre_fastp_r1="${work_dir}/three_prime_pre_fastp/mix_3/mix_3_S1_L001_R1_001.fastq.gz"
three_pre_fastp_r2="${work_dir}/three_prime_pre_fastp/mix_3/mix_3_S1_L001_R2_001.fastq.gz"
three_filtered_barcode="${work_dir}/fastp/mix3/barcode_before_reverse_complement.fastq.gz"
run_fastp_sequencer_orientation \
    "${three_pre_fastp_r1}" "${three_pre_fastp_r2}" \
    "${three_filtered_barcode}" \
    "${final_dir}/mix3_S1_L001_R2_001.fastq.gz" \
    "${work_dir}/fastp/mix3" "${CELL_LINE_FASTP_ARGS[@]}"
python3 "${REFORMAT_PY}" reverse-barcode \
    --input "${three_filtered_barcode}" \
    --output "${final_dir}/mix3_S1_L001_R1_001.fastq.gz"

python3 "${REFORMAT_PY}" split-raw \
    --raw-r1 "${CELL_LINE_RAW_R1}" \
    --raw-r2 "${CELL_LINE_RAW_R2}" \
    --indices "${SCRIPT_DIR}/configs/cell_line_indices.tsv" \
    --only-primer r \
    --outdir "${work_dir}/random_raw" --gzip-output

random_raw_r1="${work_dir}/random_raw/mix_r/mix_r_sequencer_R1.fastq.gz"
random_raw_r2="${work_dir}/random_raw/mix_r/mix_r_sequencer_R2.fastq.gz"
random_filtered_r1="${work_dir}/fastp/mixR/sequencer_R1.fastq.gz"
random_filtered_r2="${work_dir}/fastp/mixR/sequencer_R2.fastq.gz"
mkdir -p "${work_dir}/fastp/mixR"
"${fastp_command}" \
    -i "${random_raw_r1}" -I "${random_raw_r2}" \
    -o "${random_filtered_r1}" -O "${random_filtered_r2}" \
    "${CELL_LINE_RANDOM_FASTP_ARGS[@]}" \
    --json "${work_dir}/fastp/mixR/fastp.json" \
    --html "${work_dir}/fastp/mixR/fastp.html"

python3 "${REFORMAT_PY}" reformat \
    --raw-r1 "${random_filtered_r1}" \
    --raw-r2 "${random_filtered_r2}" \
    --indices "${SCRIPT_DIR}/configs/cell_line_indices.tsv" \
    --only-primer r \
    --outdir "${work_dir}/random_reformatted" \
    --tail-mode "${CELL_LINE_RANDOM_TAIL_MODE}" --gzip-output
filtered_random_r1="${work_dir}/random_reformatted/mix_r/mix_r_S1_L001_R1_001.fastq.gz"
filtered_random_r2="${work_dir}/random_reformatted/mix_r/mix_r_S1_L001_R2_001.fastq.gz"
merge_pair "${final_dir}/mixR_S1_L001_R1_001.fastq.gz" \
    "${final_dir}/mixR_S1_L001_R2_001.fastq.gz" \
    "${filtered_random_r1}" "${filtered_random_r2}"

combined_r1="${work_dir}/reformatted/mixR3_unfiltered_R1.fastq.gz"
combined_r2="${work_dir}/reformatted/mixR3_unfiltered_R2.fastq.gz"
merge_pair "${combined_r1}" "${combined_r2}" \
    "${three_r1}" "${random_r1}" "${three_r2}" "${random_r2}"
run_fastp_sequencer_orientation \
    "${combined_r1}" "${combined_r2}" \
    "${final_dir}/mixR3_S1_L001_R1_001.fastq.gz" \
    "${final_dir}/mixR3_S1_L001_R2_001.fastq.gz" \
    "${work_dir}/fastp/mixR3" "${CELL_LINE_FASTP_ARGS[@]}"

python3 "${REFORMAT_PY}" inspect "${final_dir}"/*_R1_001.fastq.gz
