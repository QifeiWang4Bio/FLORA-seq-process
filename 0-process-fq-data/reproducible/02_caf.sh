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
: "${CAF_BATCH1_RAW_R1:?Set CAF_BATCH1_RAW_R1 in the config file}"
: "${CAF_BATCH1_RAW_R2:?Set CAF_BATCH1_RAW_R2 in the config file}"
: "${CAF_BATCH2_RAW_R1:?Set CAF_BATCH2_RAW_R1 in the config file}"
: "${CAF_BATCH2_RAW_R2:?Set CAF_BATCH2_RAW_R2 in the config file}"
for path in "${CAF_BATCH1_RAW_R1}" "${CAF_BATCH1_RAW_R2}" \
    "${CAF_BATCH2_RAW_R1}" "${CAF_BATCH2_RAW_R2}"; do
    require_file "${path}"
done
require_command python3

work_dir="${FLORA_FASTQ_OUTPUT_ROOT}/caf/work"
final_dir="${FLORA_FASTQ_OUTPUT_ROOT}/caf/cellranger_fastq"
mkdir -p "${work_dir}" "${final_dir}"

for batch in batch1 batch2; do
    if [[ "${batch}" == "batch1" ]]; then
        raw_r1=${CAF_BATCH1_RAW_R1}
        raw_r2=${CAF_BATCH1_RAW_R2}
    else
        raw_r1=${CAF_BATCH2_RAW_R1}
        raw_r2=${CAF_BATCH2_RAW_R2}
    fi
    # Reproduce split_index.py followed by CutReads_3.py/Transfor_r.py.
    python3 "${REFORMAT_PY}" reformat \
        --raw-r1 "${raw_r1}" --raw-r2 "${raw_r2}" \
        --indices "${SCRIPT_DIR}/configs/caf_indices.tsv" \
        --outdir "${work_dir}/${batch}/reformatted" \
        --tail-mode trim --random-tail-mode keep \
        --gzip-output
done

for sample in sample1 sample2 sample3; do
    sample_number=S1
    for primer in 3 r; do
        library="${sample}_${primer}"
        b1="${work_dir}/batch1/reformatted/${library}/${library}_${sample_number}_L001"
        b2="${work_dir}/batch2/reformatted/${library}/${library}_${sample_number}_L001"
        merge_pair "${final_dir}/${library}_${sample_number}_L001_R1_001.fastq.gz" \
            "${final_dir}/${library}_${sample_number}_L001_R2_001.fastq.gz" \
            "${b1}_R1_001.fastq.gz" "${b2}_R1_001.fastq.gz" \
            "${b1}_R2_001.fastq.gz" "${b2}_R2_001.fastq.gz"
    done
    merge_pair "${final_dir}/${sample}_3r_${sample_number}_L001_R1_001.fastq.gz" \
        "${final_dir}/${sample}_3r_${sample_number}_L001_R2_001.fastq.gz" \
        "${final_dir}/${sample}_3_${sample_number}_L001_R1_001.fastq.gz" \
        "${final_dir}/${sample}_r_${sample_number}_L001_R1_001.fastq.gz" \
        "${final_dir}/${sample}_3_${sample_number}_L001_R2_001.fastq.gz" \
        "${final_dir}/${sample}_r_${sample_number}_L001_R2_001.fastq.gz"
done

python3 "${REFORMAT_PY}" inspect "${final_dir}"/*_R1_001.fastq.gz
