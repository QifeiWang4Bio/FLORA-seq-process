#!/usr/bin/env bash
# =============================================================================
# Module : 2-CAF / Cell Ranger count (human_sample1/2/3_3r)
# Source : breast CAF Cell Ranger notebook
# Next: 2-coverage.sh; 3-prepare-annotate.py
# =============================================================================

set -euo pipefail

sample=(sample1_3r sample2_3r sample3_3r)
for i in {0..2}
do
    echo 'Running: ' ${sample[i]}
    cd /path/to/floraseq/caf/
    /path/to/cellranger count --id=human_${sample[i]} \
    --fastqs=/path/to/floraseq/caf/Rawdata \
    --sample=${sample[i]} \
    --localmem=80 \
    --localcores=24 \
    --transcriptome=/path/to/reference/cellranger/refdata-gex-GRCh38-2020-A \
    --include-introns true
done
