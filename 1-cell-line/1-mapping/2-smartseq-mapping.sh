#!/usr/bin/env bash
# =============================================================================
# Module : 1-mapping / SMART-seq3 mapping (zUMIs)
# Dataset: HEK293T SMART-seq3 (HEK.fwdprimer, E-MTAB-8735)
#
# Final result used in this paper track (genebody / biotype / SJ diversity):
#   .../smartseq3/zumi_fwdprimer_sub10_fig/
#     Smartseq3_fwdprimer.filtered.Aligned.GeneTagged.UBcorrected.sorted.filter.bam
#
# Mapping pipeline that produced it:
#   1) download HEK.fwdprimer FASTQs
#   2) downsample to 10% (seqtk)
#   3) zUMIs -> out_dir: zumi_fwdprimer_sub10
#   4) copy/work under zumi_fwdprimer_sub10_fig for downstream figures
#
# Companion YAML: 2-smartseq3_analysis.yaml
# NOTE: absolute HPC paths — edit for your environment before running.
# =============================================================================

# conda activate sc-bio

# =============================================================================
# 0. Download HEK.fwdprimer FASTQs (E-MTAB-8735)
# =============================================================================
# Source: ftp.ebi.ac.uk/.../E-MTAB-8735/
# Only the HEK.fwdprimer library was used for the final benchmark.

cd /path/to/floraseq/benchmark/smartseq3/fq/hek_fwd

wget ftp://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment/MTAB/E-MTAB-8735/Smartseq3.HEK.fwdprimer.index1.fastq.gz .
wget ftp://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment/MTAB/E-MTAB-8735/Smartseq3.HEK.fwdprimer.index2.fastq.gz .
wget ftp://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment/MTAB/E-MTAB-8735/Smartseq3.HEK.fwdprimer.read1.fastq.gz .


# =============================================================================
# 1. Downsample to 10% (seed=100)
# =============================================================================
# Downstream figure BAM is from this subsampled run (zumi_fwdprimer_sub10*).

cd /path/to/floraseq/benchmark/smartseq3/fq/hek_fwd

seqtk sample -s 100 Smartseq3.HEK.fwdprimer.index1.fastq.gz 0.1 | pigz > sub_Smartseq3.HEK.fwdprimer.index1.fastq.gz
seqtk sample -s 100 Smartseq3.HEK.fwdprimer.index2.fastq.gz 0.1 | pigz > sub_Smartseq3.HEK.fwdprimer.index2.fastq.gz
seqtk sample -s 100 Smartseq3.HEK.fwdprimer.read1.fastq.gz 0.1 | pigz > sub_Smartseq3.HEK.fwdprimer.read1.fastq.gz


# =============================================================================
# 2. zUMIs mapping (sub10)
# =============================================================================

/path/to/software/zUMIs-main/zUMIs.sh -y 2-smartseq3_analysis.yaml


# =============================================================================
# 3. Prepare figure working directory (optional)
# =============================================================================
# Downstream genebody / biotype / SJ scripts read from zumi_fwdprimer_sub10_fig.
# After zUMIs finishes, copy (or symlink) the BAM there if needed:

# mkdir -p /path/to/floraseq/benchmark/smartseq3/zumi_fwdprimer_sub10_fig
# cp /path/to/floraseq/benchmark/smartseq3/zumi_fwdprimer_sub10/Smartseq3_fwdprimer.filtered.Aligned.GeneTagged.UBcorrected.sorted.bam \
#    /path/to/floraseq/benchmark/smartseq3/zumi_fwdprimer_sub10_fig/
