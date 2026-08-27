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
: "${MHYPO_BATCH1_RAW_R1:?Set MHYPO_BATCH1_RAW_R1 in the config file}"
: "${MHYPO_BATCH1_RAW_R2:?Set MHYPO_BATCH1_RAW_R2 in the config file}"
: "${MHYPO_BATCH2_RAW_R1:?Set MHYPO_BATCH2_RAW_R1 in the config file}"
: "${MHYPO_BATCH2_RAW_R2:?Set MHYPO_BATCH2_RAW_R2 in the config file}"
if ! declare -p MHYPO_FASTP_ARGS >/dev/null 2>&1; then
    MHYPO_FASTP_ARGS=(-q 20 -u 30 -n 10 -w 12 -l 26 --cut_right)
fi
RUN_FASTQC=${RUN_FASTQC:-1}
fastp_command=${FASTP_BIN:-fastp}
trim_galore_command=${TRIM_GALORE_BIN:-trim_galore}
fastqc_command=${FASTQC_BIN:-fastqc}

for path in "${MHYPO_BATCH1_RAW_R1}" "${MHYPO_BATCH1_RAW_R2}" \
    "${MHYPO_BATCH2_RAW_R1}" "${MHYPO_BATCH2_RAW_R2}"; do
require_file "${path}"
done
require_command python3
require_command "${fastp_command}"
require_command "${trim_galore_command}"
if [[ "${RUN_FASTQC}" == "1" ]]; then
    require_command "${fastqc_command}"
fi

work_dir="${FLORA_FASTQ_OUTPUT_ROOT}/mhypo/work"
final_dir="${FLORA_FASTQ_OUTPUT_ROOT}/mhypo/cellranger_fastq"
mkdir -p "${work_dir}" "${final_dir}"

for batch in batch1 batch2; do
    if [[ "${batch}" == "batch1" ]]; then
        raw_r1=${MHYPO_BATCH1_RAW_R1}
        raw_r2=${MHYPO_BATCH1_RAW_R2}
        three_prime_tail_mode=trim
    else
        raw_r1=${MHYPO_BATCH2_RAW_R1}
        raw_r2=${MHYPO_BATCH2_RAW_R2}
        three_prime_tail_mode=keep
    fi
    python3 "${REFORMAT_PY}" reformat \
        --raw-r1 "${raw_r1}" --raw-r2 "${raw_r2}" \
        --indices "${SCRIPT_DIR}/configs/mhypo_indices.tsv" \
        --outdir "${work_dir}/${batch}/reformatted" \
        --tail-mode keep --three-prime-tail-mode "${three_prime_tail_mode}"

    for group_number in {1..24}; do
        for primer in 3 r; do
            library="group${group_number}_${primer}"
            in_dir="${work_dir}/${batch}/reformatted/${library}"
            qc_dir="${work_dir}/${batch}/qc/${library}"
            processed_dir="${work_dir}/${batch}/processed/${library}"
            mkdir -p "${qc_dir}" "${processed_dir}"
            in_r1="${in_dir}/${library}_S1_L001_R1_001.fastq"
            in_r2="${in_dir}/${library}_S1_L001_R2_001.fastq"
            fastp_r1="${qc_dir}/${library}_S1_L001_R1_001.fastq"
            fastp_r2="${qc_dir}/${library}_S1_L001_R2_001.fastq"

            "${fastp_command}" -i "${in_r1}" -I "${in_r2}" -o "${fastp_r1}" -O "${fastp_r2}" \
                "${MHYPO_FASTP_ARGS[@]}" \
                --json "${qc_dir}/fastp.json" --html "${qc_dir}/fastp.html"

            final_qc_r1=${fastp_r1}
            final_qc_r2=${fastp_r2}
            if [[ "${primer}" == "r" ]]; then
                "${trim_galore_command}" -q 20 --phred33 --stringency 3 --length 26 -e 0.1 \
                    --paired --no_report_file "${fastp_r1}" "${fastp_r2}" -o "${qc_dir}"
                pass1_r1="${qc_dir}/${library}_S1_L001_R1_001_val_1.fq"
                pass1_r2="${qc_dir}/${library}_S1_L001_R2_001_val_2.fq"
                "${trim_galore_command}" -q 20 --phred33 --illumina --stringency 3 --length 26 -e 0.1 \
                    --paired --no_report_file "${pass1_r1}" "${pass1_r2}" -o "${qc_dir}"
                final_qc_r1="${qc_dir}/${library}_S1_L001_R1_001_val_1_val_1.fq"
                final_qc_r2="${qc_dir}/${library}_S1_L001_R2_001_val_2_val_2.fq"
            fi

            if [[ "${RUN_FASTQC}" == "1" ]]; then
                "${fastqc_command}" -t 12 -o "${qc_dir}" "${final_qc_r1}" "${final_qc_r2}"
            fi

            out_r1="${processed_dir}/${library}_S1_L001_R1_001.fastq.gz"
            out_r2="${processed_dir}/${library}_S1_L001_R2_001.fastq.gz"
            python3 "${REFORMAT_PY}" truncate --input "${final_qc_r1}" \
                --output "${out_r1}" --length 26
            temp_r1="${processed_dir}/${library}_temporary_R1.fastq.gz"
            merge_pair "${temp_r1}" "${out_r2}" "${out_r1}" "${final_qc_r2}"
            mv "${temp_r1}" "${out_r1}"
        done
    done
done

for group_number in {1..24}; do
    for primer in 3 r; do
        library="group${group_number}_${primer}"
        b1="${work_dir}/batch1/processed/${library}/${library}_S1_L001"
        b2="${work_dir}/batch2/processed/${library}/${library}_S1_L001"
        merge_pair "${final_dir}/${library}_S1_L001_R1_001.fastq.gz" \
            "${final_dir}/${library}_S1_L001_R2_001.fastq.gz" \
            "${b1}_R1_001.fastq.gz" "${b2}_R1_001.fastq.gz" \
            "${b1}_R2_001.fastq.gz" "${b2}_R2_001.fastq.gz"
    done
    merge_pair "${final_dir}/group${group_number}_3r_S1_L001_R1_001.fastq.gz" \
        "${final_dir}/group${group_number}_3r_S1_L001_R2_001.fastq.gz" \
        "${final_dir}/group${group_number}_3_S1_L001_R1_001.fastq.gz" \
        "${final_dir}/group${group_number}_r_S1_L001_R1_001.fastq.gz" \
        "${final_dir}/group${group_number}_3_S1_L001_R2_001.fastq.gz" \
        "${final_dir}/group${group_number}_r_S1_L001_R2_001.fastq.gz"
done

python3 "${REFORMAT_PY}" inspect "${final_dir}"/*_R1_001.fastq.gz
