#!/usr/bin/env bash
# =============================================================================

# -----------------------------------------------------------------------------
# Seed content from 1-benchmark-data.sh
# -----------------------------------------------------------------------------

## 1. FLORA-seq
# Library designs:
#   mixR3 = 3'+Random (3+R)
#   mix3  = 3' only
#   mixR  = Random only
# Reference modes:
#   dual   = hg38 + mm39 (species-mix)
#   human  = GRCh38 only
#   mouse  = GRCm39 only
cd /path/to/floraseq/cell_line

# ---------------------------------------------------------------------------
# 1.1 Dual-genome mapping (hg38 + mm39)
# ---------------------------------------------------------------------------

# 3'+Random (mixR3), dual
cellranger count --id=mixR3 \
--fastqs=/path/to/floraseq/cell_line/Cleandata \
--sample=mixR3 \
--localmem=80 \
--localcores=24 \
--transcriptome=/path/to/reference/cellranger/hg38_and_mm39 \
--include-introns true

# 3' only (mix3), dual
cellranger count --id=mix3 \
--fastqs=/path/to/floraseq/cell_line/Cleandata \
--sample=mix3 \
--localmem=80 \
--localcores=16 \
--transcriptome=/path/to/reference/cellranger/hg38_and_mm39 \
--include-introns true

# Random only (mixR), dual
cellranger count --id=mixR \
--fastqs=/path/to/floraseq/cell_line/Cleandata \
--sample=mixR \
--localmem=80 \
--localcores=24 \
--transcriptome=/path/to/reference/cellranger/hg38_and_mm39 \
--include-introns true

# ---------------------------------------------------------------------------
# 1.2 Human-only mapping (GRCh38)
# ---------------------------------------------------------------------------

# 3'+Random (mixR3), human
cellranger count --id=mixR3_human \
--fastqs=/path/to/floraseq/cell_line/Cleandata \
--sample=mixR3 \
--localmem=200 \
--localcores=24 \
--transcriptome=/path/to/reference/cellranger/refdata-gex-GRCh38-2020-A \
--include-introns true

# 3' only (mix3), human
cellranger count --id=mix3_human \
--fastqs=/path/to/floraseq/cell_line/Cleandata \
--sample=mix3 \
--localmem=80 \
--localcores=16 \
--transcriptome=/path/to/reference/cellranger/refdata-gex-GRCh38-2020-A \
--include-introns true

# Random only (mixR), human
cellranger count --id=mixR_human \
--fastqs=/path/to/floraseq/cell_line/Cleandata \
--sample=mixR \
--localmem=80 \
--localcores=16 \
--transcriptome=/path/to/reference/cellranger/refdata-gex-GRCh38-2020-A \
--include-introns true

# ---------------------------------------------------------------------------
# 1.3 Mouse-only mapping (GRCm39)
# ---------------------------------------------------------------------------

# 3' only (mix3), mouse
cellranger count --id=mix3_mouse \
--fastqs=/path/to/floraseq/cell_line/Cleandata \
--sample=mix3 \
--localmem=80 \
--localcores=16 \
--transcriptome=/path/to/reference/cellranger/GRCm39 \
--include-introns true

# Random only (mixR), mouse
cellranger count --id=mixR_mouse \
--fastqs=/path/to/floraseq/cell_line/Cleandata \
--sample=mixR \
--localmem=80 \
--localcores=16 \
--transcriptome=/path/to/reference/cellranger/GRCm39 \
--include-introns true



# ===================================================================================================================
## 2. 10x-3': 293T + NIH3T3 (1k cell) - 3' v3.1:
# https://www.10xgenomics.com/datasets/1-k-1-1-mixture-of-human-hek-293-t-and-mouse-nih-3-t-3-cells-3-v-3-1-3-1-standard-6-0-0
cd /path/to/floraseq/benchmark/10x_3/293t_mix_1k_v3-1

# human reads
echo "Job started on $(date)"
cellranger count \
  --id=293t_mix_1k_v3-1_human \
  --sample=1k_hgmm_3p \
  --fastqs=/path/to/floraseq/benchmark/10x_3/293t_mix_1k_v3-1/1k_hgmm_3p_fastqs \
  --transcriptome=/path/to/reference/cellranger/refdata-gex-GRCh38-2020-A \
  --include-introns=true \
  --localcores=16 \
  --localmem=150

echo "Job finished on $(date)"

# human+mouse reads
echo "Job started on $(date)"
cellranger count \
  --id=293t_mix_1k_v3-1_run \
  --sample=1k_hgmm_3p \
  --fastqs=/path/to/floraseq/benchmark/10x_3/293t_mix_1k_v3-1/1k_hgmm_3p_fastqs \
  --transcriptome=/path/to/reference/cellranger/hg38_and_mm39 \
  --include-introns=true \
  --localcores=16 \
  --localmem=150

echo "Job finished on $(date)"


# ===================================================================================================================
## 10x 5'
/path/to/floraseq/benchmark/bin/cellranger_10x_5.sh
#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=150G
#SBATCH --time=48:00:00
#SBATCH --partition=corexd192
#SBATCH --job-name=wang_
#SBATCH --output=log/cellranger_mix_10x_5_human.out
#SBATCH --error=log/cellranger_mix_10x_5_human.err
cd /path/to/floraseq/benchmark/10x_5

# human
echo "Job started on $(date)"
cellranger count \
  --id=10k_hgmm_5pv2_human \
  --sample=10k_hgmm_5pv2_nextgem_Chromium_X \
  --fastqs=/path/to/floraseq/benchmark/10x_5/10k_hgmm_5pv2_nextgem_Chromium_X_fastqs \
  --transcriptome=/path/to/reference/cellranger/refdata-gex-GRCh38-2020-A \
  --include-introns=true \
  --localcores=16 \
  --localmem=150

echo "Job finished on $(date)"

# hg38+mm39
echo "Job started on $(date)"
cellranger count \
  --id=10k_hgmm_5pv2_run \
  --sample=10k_hgmm_5pv2_nextgem_Chromium_X \
  --fastqs=/path/to/floraseq/benchmark/10x_5/10k_hgmm_5pv2_nextgem_Chromium_X_fastqs \
  --transcriptome=/path/to/reference/cellranger/hg38_and_mm39 \
  --include-introns=true \
  --localcores=16 \
  --localmem=150

echo "Job finished on $(date)"


# ===================================================================================================================
## VASA-seq
cd /path/to/floraseq/benchmark/vasaseq/fq

# Rename Read 1 (usually the Barcode/UMI read)
mv SRR14783059_1.fastq.gz SRR14783059_S1_L001_R1_001.fastq.gz
# Rename Read 2 (usually the cDNA read)
mv SRR14783059_2.fastq.gz SRR14783059_S1_L001_R2_001.fastq.gz

cd /path/to/floraseq/benchmark/vasaseq

# Remove the 3' poly(A) tail from Read 2 with cutadapt.
# -a "A{10}" trims a terminal run of ten A bases.
# -m 20 discards reads shorter than 20 bp after trimming.
# -j 8 uses eight threads.

cutadapt -a "A{10}" -m 20 -j 16 \
-o SRR14783059_R2_trimmed.fastq.gz \
SRR14783059_S1_L001_R2_001.fastq.gz

# Build a new STAR genome index.
cd /path/to/reference/star_index_2.7.11b
STAR --runMode genomeGenerate \
     --runThreadN 16 \
     --genomeDir /path/to/reference/star_index_2.7.11b \
     --genomeFastaFiles /path/to/reference/cellranger/refdata-gex-GRCh38-2020-A/fasta/genome.fa \
     --sjdbGTFfile /path/to/reference/cellranger/refdata-gex-GRCh38-2020-A/genes/genes.gtf \
     --sjdbOverhang 100

# Example STAR index: /path/to/STAR_Index
# Example GTF: /path/to/genes.gtf (recommended when building the index).

/path/to/floraseq/benchmark/bin/star_vasaseq.sh
#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=150G
#SBATCH --time=48:00:00
#SBATCH --partition=corexd192
#SBATCH --job-name=wang_
#SBATCH --output=log/star_vasaseq.out
#SBATCH --error=log/star_vasaseq.err
source /path/to/software/miniforge3/etc/profile.d/conda.sh
conda activate sc-bio

cd /path/to/floraseq/benchmark/vasaseq
echo "Job started on $(date)"
# Retain reads mapped to intronic regions in the BAM file.
STAR --runThreadN 16 \
     --genomeDir /path/to/reference/star_2.7 \
     --readFilesIn /path/to/floraseq/benchmark/vasaseq/fq/SRR14783059_R2_trimmed.fastq.gz \
     --readFilesCommand zcat \
     --outFileNamePrefix vasa_293t \
     --outSAMtype BAM SortedByCoordinate \
     --outSAMunmapped Within \
     --outSAMattributes NH HI AS nM NM MD \
     --limitBAMsortRAM 30000000000

echo "Job finished on $(date)"


# ===================================================================================================================
## SMART-seq 
# test
cd /path/to/floraseq/benchmark/smartseq3/hek293t_coverage
#mkdir -p /path/to/floraseq/benchmark/smartseq3/hek293t_coverage/star_tmp_smartseq3

/path/to/software/miniforge3/envs/sc-py/bin/STAR \
--runThreadN 24 \
--runMode genomeGenerate \
--genomeDir /path/to/reference/star_index_2.7.4 \
--genomeFastaFiles /path/to/reference/cellranger/refdata-gex-GRCh38-2020-A/fasta/genome.fa \
--sjdbGTFfile /path/to/reference/cellranger/refdata-gex-GRCh38-2020-A/genes/genes.gtf \
--sjdbOverhang 74

STAR \
  --runThreadN 16 \
  --genomeDir /path/to/reference/star_index_2.7.4 \
  --readFilesIn \
    /path/to/floraseq/benchmark/smartseq3/fq/hek_fwd/sub_Smartseq3.HEK.fwdprimer.read1.fastq.gz \
  --readFilesCommand zcat \
  --outFileNamePrefix /path/to/floraseq/benchmark/smartseq3/hek293t_coverage/star \
  --outSAMtype BAM SortedByCoordinate \
  --outSAMunmapped Within \
  --outSAMattributes NH HI AS nM NM MD \
  --limitBAMsortRAM 30000000000 \
  --outTmpDir /path/to/floraseq/benchmark/smartseq3/hek293t_coverage/star_tmp_smartseq3


