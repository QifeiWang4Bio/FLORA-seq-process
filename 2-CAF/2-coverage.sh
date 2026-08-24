#!/usr/bin/env bash
# =============================================================================
# Module : 2-CAF / BAM filter + coverage / reads track (step 1)
# Notebook-style sections — run selectively after configuring input paths.
# Prerequisite: 1-process-caf.sh (or existing Cell Ranger outs)
# =============================================================================

## ======================================================================
# 1. filter umapped reads
# 3'
sample=(human_sample1_3 human_sample2_3 human_sample3_3)
for i in {0..2}
do
    echo 'Running: ' ${sample[i]}
    cd /path/to/floraseq/caf/${sample[i]}/outs
    samtools view -b \
    -F 4 -F 256 -F 2048 \
    possorted_genome_bam.bam > possorted_genome_bam_filter.bam
    samtools index possorted_genome_bam_filter.bam
done


# Random
#!/bin/bash
sample=(human_sample1_r human_sample2_r human_sample3_r)
for i in {0..2}
do
    echo 'Running: ' ${sample[i]}
    cd /path/to/floraseq/caf/${sample[i]}/outs
    samtools view -b \
    -F 4 -F 256 -F 2048 \
    possorted_genome_bam.bam > possorted_genome_bam_filter.bam
    samtools index possorted_genome_bam_filter.bam
done

# ================================================================================================
# 2. reads balance
# 3' -> 10%
#!/bin/bash
sample=(human_sample1_3 human_sample2_3 human_sample3_3)
for i in {0..2}
do
    echo 'Running: ' ${sample[i]}
    cd /path/to/floraseq/caf/${sample[i]}/outs
    samtools view -b -s 0.1 possorted_genome_bam_filter.bam > possorted_genome_bam_filter_${sample[i]}.bam
    samtools index possorted_genome_bam_filter_${sample[i]}.bam
done

# Random -> 20%
#!/bin/bash
sample=(human_sample1_r human_sample2_r human_sample3_r)
for i in {0..2}
do
    echo 'Running: ' ${sample[i]}
    cd /path/to/floraseq/caf/${sample[i]}/outs
    samtools view -b -s 0.2 possorted_genome_bam_filter.bam > possorted_genome_bam_filter_${sample[i]}.bam
    samtools index possorted_genome_bam_filter_${sample[i]}.bam
done


# ================================================================================================
# 3. reads coverage
: <<'BAM_PATHS'
/path/to/floraseq/caf/human_sample1_3/outs/possorted_genome_bam_filter_human_sample1_3.bam
/path/to/floraseq/caf/human_sample2_3/outs/possorted_genome_bam_filter_human_sample2_3.bam
/path/to/floraseq/caf/human_sample3_3/outs/possorted_genome_bam_filter_human_sample3_3.bam
/path/to/floraseq/caf/human_sample1_r/outs/possorted_genome_bam_filter_human_sample1_r.bam
/path/to/floraseq/caf/human_sample2_r/outs/possorted_genome_bam_filter_human_sample2_r.bam
/path/to/floraseq/caf/human_sample3_r/outs/possorted_genome_bam_filter_human_sample3_r.bam
BAM_PATHS

# /path/to/floraseq/caf/bin/cover_caf.sh
#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=50G
#SBATCH --time=99:00:00
#SBATCH --partition=corexd192
#SBATCH --job-name=wang_
#SBATCH --output=log/genebody_caf.out
#SBATCH --error=log/genebody_caf.err
source /path/to/software/miniforge3/etc/profile.d/conda.sh
conda activate sc-bio
echo "Job started on $(date)"

cd /path/to/floraseq/caf/coverage
/path/to/software/miniforge3/envs/sc-bio/bin/geneBody_coverage.py \
-r /path/to/reference/RSeQC_ref/hg38_GENCODE.v38.bed \
-i cover_caf.txt \
-o RseQC_caf

echo "Job finished on $(date)"



# ================================================================================================
# 4. reads coverage - fig by python


# Output: RseQC_caf.geneBodyCoverage.txt

python3 - <<'PYTHON'
import matplotlib.pyplot as plt
import pandas as pd

cover_caf = pd.read_csv('/path/to/floraseq/caf/coverage/RseQC_caf.geneBodyCoverage.txt',sep='\t',index_col=0)

cover_caf.index = ['sample1-3', 'sample2-3', 'sample3-3', 'sample1-r', 'sample2-r', 'sample3-r']
cover_caf.to_csv('/path/to/floraseq/caf/coverage/caf_genebody_coverage.csv', sep='\t')

# Calculate proportions
cover_prop = cover_caf.div(cover_caf.sum(axis=1), axis=0)
plt.rcParams['pdf.fonttype']=42


# Plot gene body coverage
plt.figure(figsize=(10, 6))
plt.grid(False) 
for tech in cover_prop.index:
    alpha_val = 1.0 if tech in ['FLORA-seq-3', 'FLORA-seq-R'] else 0.8
    plt.plot(cover_prop.columns, cover_prop.loc[tech], label=tech, linewidth=2, alpha=alpha_val)

plt.xlabel('Gene body percentile (%)', fontsize=12)
plt.ylabel('Coverage proportion', fontsize=12)
plt.title('Gene Body Coverage Distribution', fontsize=14)
plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
plt.xticks([1, 11, 21, 31, 41, 51, 61, 71, 81, 91, 100], [0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100])
plt.tight_layout()
plt.savefig('/path/to/floraseq/caf/coverage/caf_genebody_coverage.pdf', dpi=300, bbox_inches='tight')
plt.close()
PYTHON




