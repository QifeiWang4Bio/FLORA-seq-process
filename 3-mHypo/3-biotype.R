#!/usr/bin/env Rscript
# =============================================================================
# Module : 3-mHypo / Gene biotype fractions (RCAS) on merge_3r.bed
# Source : 3_mHypo/bin/old/biotype_lncRNA.R
# Prerequisite: merge_3r.bed (BAM→BED not in this script; see KNOWN_ISSUES.md)
# Environment: R with RCAS, data.table, ggplot2
# =============================================================================

# bam to bed

library(RCAS)
gff <- importGtf(filePath = '/path/to/reference/ensembl/mouse/Mus_musculus.GRCm39.108.gtf', keepStandardChr = FALSE)
# lnc_pos = which(gff$gene_biotype == 'lncRNA')
# lnc_gene = gff$gene_name[lnc_pos]
# write.table(lnc_gene, '/path/to/reference/ensembl/mouse/lncRNA_gene.txt',quote =FALSE,col.names=FALSE, row.names=FALSE)

setwd('/path/to/floraseq/mhypo/merge_3r_')
file = 'merge_3r.bed'
# 
queryRegions <- importBed(filePath = file, keepStandardChr = FALSE)
queryRegions <- sample(queryRegions, length(queryRegions) * 0.2)   # Downsample because the query count exceeds the valid queryIndex count.


# 
# 1. Querying the annotation file
overlaps <- as.data.table(queryGff(queryRegions = queryRegions, gffData = gff))   # 599844412
overlaps_multi <- overlaps[,c(14,28)]     # 14: biotype   28: queryIndex
overlaps_multi_unique <- unique(overlaps_multi)   # A query index may still map to multiple biotypes.
freq <- as.data.frame(table(overlaps_multi_unique$queryIndex))  # Values above one identify indices assigned to multiple biotypes.
multi_index <- as.numeric(as.character(freq[freq$Freq>1, ]$Var1))   # Extract ambiguous indices.
print(paste0('Multi annotated reads: ', length(multi_index)))   # 1659372
# 
# 2. Finding targeted gene types
biotype_col <- grep('gene_biotype', colnames(overlaps), value = T)
overlaps_uniq <- overlaps[!duplicated(overlaps$queryIndex), ]   # Keep one provisional biotype per read index.
rm_index <- match(multi_index, overlaps_uniq$queryIndex)    # Locate one-to-many indices for removal.
overlaps_uniq_rmMulti <- overlaps_uniq[-rm_index, ]
print(paste0('Overlap to gtf: ', length(unique(overlaps_uniq$queryIndex))))  # ! 75541680
# 
df <- overlaps_uniq_rmMulti[, length(unique(queryIndex)), by=biotype_col]
df <- data.frame(df)
rownames(df) <- df$gene_biotype
write.table(df, '/path/to/floraseq/mhypo/seurat_out/mHypo_meng_hvg5k_C25/lncRNA/Biotype_overlap.txt',sep=',', row.names=FALSE)

# plot by Biotype_overlap.tx

