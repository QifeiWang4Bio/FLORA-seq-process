#!/usr/bin/env bash

set -euo pipefail

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
REFORMAT_PY="${SCRIPT_DIR}/flora_fastq_reformat.py"

require_file() {
    local path=$1
    if [[ ! -f "${path}" ]]; then
        echo "Required input file not found: ${path}" >&2
        exit 1
    fi
}

require_command() {
    local command_name=$1
    if ! command -v "${command_name}" >/dev/null 2>&1; then
        echo "Required command not found: ${command_name}" >&2
        exit 1
    fi
}

run_fastp_sequencer_orientation() {
    local barcode_r1=$1
    local cdna_r2=$2
    local output_r1=$3
    local output_r2=$4
    local report_dir=$5
    shift 5
    local fastp_command=${FASTP_BIN:-fastp}
    local temporary_cdna="${report_dir}/fastp_cDNA_R1.fastq.gz"
    local temporary_barcode="${report_dir}/fastp_barcode_R2.fastq.gz"

    mkdir -p "${report_dir}" "$(dirname "${output_r1}")"
    "${fastp_command}" \
        -i "${cdna_r2}" -I "${barcode_r1}" \
        -o "${temporary_cdna}" -O "${temporary_barcode}" \
        "$@" \
        --json "${report_dir}/fastp.json" \
        --html "${report_dir}/fastp.html"
    mv "${temporary_barcode}" "${output_r1}"
    mv "${temporary_cdna}" "${output_r2}"
}

merge_pair() {
    local output_r1=$1
    local output_r2=$2
    shift 2
    local inputs=("$@")
    local half=$((${#inputs[@]} / 2))
    python3 "${REFORMAT_PY}" merge \
        --r1 "${inputs[@]:0:${half}}" \
        --r2 "${inputs[@]:${half}}" \
        --output-r1 "${output_r1}" \
        --output-r2 "${output_r2}"
}
