#!/usr/bin/env Rscript
# =============================================================================
# Module : 2-CAF / MARVEL SJ analysis (CAF subtypes)
# Dataset: Breast cancer FLORA-seq CAF cells (3 samples, 3'+Random)
# Source : 2_hBreast/bin/5-sj-marvel.R
#
# Prerequisite:
#   1) 7-sj-process.sh  → STARsolo Gene/SJ under caf_sc_update/human_sample*_3r
#   2) 6-caf-analysis.py → meta.csv (cell_id, sample, barcode, celltype)
#
# Final path used in the paper track:
#   - Merge human_sample1/2/3_3r
#   - Compare mCAF vs tCAF
#   - Export DE tables + gene-level SJ figures
#   - Save marvel_CAF_3r.rds
#
# Later sections (candidate gene filters / TPM2 / length filters) are
# exploratory follow-ups; keep or skip as needed.
#
# Sibling scripts NOT used here:
#   5-sj-marvel-merge.R   — same merge + downsample; mCAF vs dCAF
#   5-sj-marvel-single.R  — single-sample trial
# =============================================================================

# library(tidyverse)
library(Seurat)
library(R.utils)
library(Matrix)
library(magrittr)
library(ggplot2)

WORK_DIR <- "/path/to/floraseq/marvel/caf_sc_update"
META_CSV <- "/path/to/floraseq/marvel/CAF_sc/meta.csv"
GTF_FILE <- "/path/to/floraseq/marvel/basic/gtf/gencode.v31.annotation.gtf"

setwd(WORK_DIR)
meta <- read.table(META_CSV, sep = "\t", header = TRUE, row.names = 1)

# Remap meta rownames when analyzing _3 or _r libraries (no-op for 3r).
update_meta_prefix <- function(meta, sample_suffix) {
  if (sample_suffix == "3r") {
    return(meta)
  }
  # ^(human_sample\d+)_([^_]+)_  → keep sample id, replace library suffix
  pattern <- "^(human_sample\\d+)_([^_]+)_"
  replacement <- paste0("\\1_", sample_suffix, "_")
  rownames(meta) <- sub(pattern = pattern, replacement = replacement, x = rownames(meta))
  meta$cell_id.1 <- sub(pattern = pattern, replacement = replacement, x = meta$cell_id.1)
  return(meta)
}

# FINAL setting: three samples, 3'+Random
sample_suffix <- "3r"
samples <- paste0("human_sample", seq(1, 3), "_", sample_suffix)
output_base <- paste0("./sample_", sample_suffix)

meta <- update_meta_prefix(meta, sample_suffix = sample_suffix)

# =============================================================================
# 1. Process gene expression matrix (STARsolo Gene → MARVEL input)
# =============================================================================
dir.create(output_base, showWarnings=FALSE, recursive=TRUE)
dir.create(paste0(output_base, "/gene/raw"), showWarnings=FALSE, recursive=TRUE)
dir.create(paste0(output_base, "/gene/norm"), showWarnings=FALSE, recursive=TRUE)

# Raw gene count
sce.list <- lapply(samples, function(x){
    cat(paste0("\nProcessing: ", x))
    fls <- list.files(paste0(x, "/Solo.out/Gene/filtered"), full.names=TRUE)
    if(length(grep("gz", fls)) != 3){
        for(fl in fls){ gzip(fl, remove=FALSE) }
    }
    count <- Read10X(paste0(x, "/Solo.out/Gene/filtered"))
    sce <- CreateSeuratObject(count)
    sce <- RenameCells(sce, add.cell.id=x)
    sce$Group <- x
    bar <- intersect(colnames(sce), rownames(meta))
    sce <- sce[, bar]
    print(paste("Cell counts:", dim(sce)[2]))
    return(sce)
    })
names(sce.list) <- samples

# Merge samples
sce <- merge(sce.list[[1]], sce.list[-1])
table(sce$Group)
count_merge <- sce@assays$RNA@counts
dim(count_merge)

writeMM(count_merge, file=paste0(output_base, "/gene/raw/matrix.mtx"))


# Barcode
count_merge_phe <- sce@meta.data %>% 
    tibble::rownames_to_column("cell.id") %>% 
    dplyr::select(cell.id)
write.table(count_merge_phe, sep="\t", row.names=FALSE, 
    file=paste0(output_base, "/gene/raw/phenoData.txt"))

# Gene features
count_merge_feat <- data.frame(gene_short_name=rownames(sce))
write.table(count_merge_feat, sep="\t", row.names=FALSE, 
    file=paste0(output_base, "/gene/raw/featureData.txt"))

# Normalized gene count
sce <- NormalizeData(sce)

count_merge_norm <- Matrix::Matrix(exp(sce@assays$RNA@data)-1, sparse=TRUE)
count_merge_norm[1:4, 1:4]
writeMM(count_merge_norm, file=paste0(output_base, "/gene/norm/matrix_normalised.mtx"))

# Barcode with group info
count_merge_norm_phe <- sce@meta.data %>% 
    tibble::rownames_to_column("cell.id") %>% 
    dplyr::select(cell.id, Group)
write.table(count_merge_norm_phe, sep="\t", row.names=FALSE, 
    file=paste0(output_base, "/gene/norm/phenoData.txt"))

# Gene features
count_merge_norm_feat <- data.frame(gene_short_name=rownames(sce))
write.table(count_merge_norm_feat, sep="\t", row.names=FALSE, 
    file=paste0(output_base, "/gene/norm/featureData.txt"))


# =============================================================================
# 2. Process SJ matrix (STARsolo SJ → MARVEL input)
# =============================================================================
dir.create(paste0(output_base, "/sj/raw"), showWarnings=FALSE, recursive=TRUE)

library(data.table)
library(dplyr)

# Build features.tsv from SJ.out.tab for each sample
for (x in samples) {
  print(paste("Processing sample:", x))
  
  sj_file <- paste0(x, "/SJ.out.tab")
  print(paste("Looking for:", sj_file))
  
  if (!file.exists(sj_file)) {
    warning("SJ.out.tab not found for ", x)
    print(paste("File does not exist:", sj_file))
    next
  }
  
  print(paste("Found SJ.out.tab for", x))
  
  sj_dir <- paste0(x, "/Solo.out/SJ/raw")
  dir.create(sj_dir, showWarnings = FALSE, recursive = TRUE)
  print(paste("Created directory:", sj_dir))
  
  sj_name <- data.table::fread(sj_file, header = FALSE, data.table = FALSE) %>% 
    dplyr::mutate(coord.intron = paste(V1, V2, V3, sep = ":")) %>% 
    dplyr::select(coord.intron)
  
  print(paste("Read", nrow(sj_name), "SJ coordinates from SJ.out.tab"))
  print(paste("Example coordinates:", head(sj_name$coord.intron, 3)))
  
  output_file <- paste0(sj_dir, "/features.tsv")
  if(file.exists(output_file)){
    file.remove(output_file)
    print(paste("Removed existing file:", output_file))
  }
  
  write.table(
    sj_name,
    sep = "\t",
    row.names = FALSE,
    col.names = FALSE,
    quote = FALSE,
    file = output_file
  )
  print(paste("Written features.tsv to:", output_file))
  print(paste("File exists:", file.exists(output_file)))
}


# Load SJ counts (long time)
sce_sj.list <- lapply(samples, function(x){
    print(x)
    fls <- list.files(paste0(x, "/Solo.out/SJ/raw"), full.names=TRUE)
    if(length(grep("gz", fls)) != 3){
        for(fl in fls){ gzip(fl, remove=FALSE) }
    }
    count <- Read10X(paste0(x, "/Solo.out/SJ/raw"), gene.column=1)
    if(is.list(count)){ count <- count[[1]] }
    sce_sj <- CreateSeuratObject(count)
    sce_sj <- RenameCells(sce_sj, add.cell.id=x)
    sce_sj <- sce_sj[, colnames(sce.list[[x]])]
    return(sce_sj)
})

sce_sj <- merge(sce_sj.list[[1]], sce_sj.list[-1])
sj_count_merge <- sce_sj@assays$RNA@counts
writeMM(sj_count_merge, file=paste0(output_base, "/sj/raw/matrix.mtx"))

sj_count_merge_phe <- sce_sj@meta.data %>% 
    tibble::rownames_to_column("cell.id") %>% 
    dplyr::select(cell.id)
write.table(sj_count_merge_phe, sep="\t", row.names=FALSE, 
    file=paste0(output_base, "/sj/raw/phenoData.txt"))

sj_count_merge_feat <- data.frame(coord.intron=rownames(sce_sj)) %>% 
    dplyr::mutate(coord.intron=gsub("-", "_", coord.intron))
write.table(sj_count_merge_feat, sep="\t", row.names=FALSE, 
    file=paste0(output_base, "/sj/raw/featureData.txt"))


# =============================================================================
# 3. Cell dimension reduction (UMAP for MARVEL)
# =============================================================================
dir.create(paste0(output_base, "/umap"), showWarnings=FALSE, recursive=TRUE)

sce <- sce %>%
    NormalizeData() %>%
    FindVariableFeatures() %>%
    ScaleData() %>% 
    Seurat::RunPCA() %>% 
    Seurat::RunUMAP(dims=1:30)

umap_embed <- sce@reductions$umap@cell.embeddings %>% 
    as.data.frame() %>% 
    tibble::rownames_to_column("cell.id") %>%
    dplyr::rename(x=UMAP_1, y=UMAP_2) 
write.table(umap_embed, row.names=FALSE, sep="\t",
    file=paste0(output_base, "/umap/UMAP_coordinate.txt"))


# =============================================================================
# 4. Create MARVEL object
# =============================================================================
library(MARVEL)

setwd(WORK_DIR)

df.gene.norm <- readMM(paste0(output_base, "/gene/norm/matrix_normalised.mtx"))
df.gene.norm.pheno <- read.table(paste0(output_base, "/gene/norm/phenoData.txt"), sep="\t", header=TRUE)
df.gene.norm.feature <- read.table(paste0(output_base, "/gene/norm/featureData.txt"), sep="\t", header=TRUE)

df.gene.count <- readMM(paste0(output_base, "/gene/raw/matrix.mtx"))
df.gene.count.pheno <- read.table(paste0(output_base, "/gene/raw/phenoData.txt"), sep="\t", header=TRUE)
df.gene.count.feature <- read.table(paste0(output_base, "/gene/raw/featureData.txt"), sep="\t", header=TRUE)

df.sj.count <- readMM(paste0(output_base, "/sj/raw/matrix.mtx"))
df.sj.count.pheno <- read.table(paste0(output_base, "/sj/raw/phenoData.txt"), sep="\t", header=TRUE)
df.sj.count.feature <- read.table(paste0(output_base, "/sj/raw/featureData.txt"), sep="\t", header=TRUE)

row.names(df.sj.count) <- df.sj.count.feature$coord.intron
df.sj.count.feature$coord.intron <- row.names(df.sj.count)
rownames(df.sj.count) <- df.sj.count.feature[, 1]
colnames(df.sj.count) <- df.sj.count.pheno[, 1]

df.coord <- read.table(paste0(output_base, "/umap/UMAP_coordinate.txt"), sep="\t", header=TRUE)
meta <- meta[colnames(df.sj.count), ]

gtf <- data.table::fread(GTF_FILE, sep = "\t", header = FALSE) %>% as.data.frame()
gtf$V9 <- gsub("gene_type", "gene_biotype", gtf$V9)

# Create MARVEL object
marvel <- CreateMarvelObject.10x(
    gene.norm.matrix=df.gene.norm,
    gene.norm.pheno=df.gene.norm.pheno,
    gene.norm.feature=df.gene.norm.feature,
    gene.count.matrix=df.gene.count,
    gene.count.pheno=df.gene.count.pheno,
    gene.count.feature=df.gene.count.feature,
    sj.count.matrix=df.sj.count,
    sj.count.pheno=df.sj.count.pheno,
    sj.count.feature=df.sj.count.feature,
    pca=df.coord,
    gtf=gtf
)


# =============================================================================
# 5. Annotate / validate SJ; attach CAF subtypes from meta.csv
# =============================================================================
marvel <- AnnotateGenes.10x(MarvelObject=marvel)
head(marvel$gene.metadata)
table(marvel$gene.metadata$gene_type)

marvel <- AnnotateSJ.10x(MarvelObject=marvel)
head(marvel$sj.metadata)

marvel <- ValidateSJ.10x(MarvelObject=marvel)
head(marvel$sj.metadata)

marvel <- FilterGenes.10x(MarvelObject=marvel,
                          gene.type="protein_coding"
                          )
head(marvel$sj.metadata)

# Add cell type information
rownames(marvel$sample.metadata) <- marvel$sample.metadata$cell.id
marvel$sample.metadata$celltype <- meta$celltype

# Create cell group list
cell.group.list <- split(marvel$sample.metadata$cell.id, 
                        marvel$sample.metadata$celltype)


# =============================================================================
# 6. Expression / splicing percentage (FINAL contrast: mCAF vs tCAF)
# =============================================================================
celltype_1 = "mCAF"
celltype_2 = "tCAF"
celltype_dir <- paste0(celltype_1, "_", celltype_2)
output_celltype <- paste0(output_base, "/", celltype_dir)
dir.create(output_celltype, showWarnings=FALSE, recursive=TRUE)

marvel <- PlotPctExprCells.Genes.10x(
    MarvelObject=marvel,
    cell.group.g1=cell.group.list[[celltype_1]],
    cell.group.g2=cell.group.list[[celltype_2]],
    min.pct.cells=5
)
head(marvel$pct.cells.expr$Gene$Data)
p1 <- marvel$pct.cells.expr$Gene$Plot
ggsave(p1, filename=paste0(output_celltype, "/fig1_exp.png"))

# SJ percentage
marvel <- PlotPctExprCells.SJ.10x(
    MarvelObject=marvel,
    cell.group.g1=cell.group.list[[celltype_1]],
    cell.group.g2=cell.group.list[[celltype_2]],
    min.pct.cells.genes=5,
    min.pct.cells.sj=5,
    downsample=TRUE,
    downsample.pct.sj=100
)
head(marvel$pct.cells.expr$SJ$Data)
p2 <- marvel$pct.cells.expr$SJ$Plot
ggsave(p2, filename=paste0(output_celltype, "/fig2_sj.png"))


# =============================================================================
# 7. Differential SJ analysis + isoform switch
# =============================================================================

# Compare SJ between groups
marvel <- CompareValues.SJ.10x(
    MarvelObject=marvel,
    cell.group.g1=cell.group.list[[celltype_1]],
    cell.group.g2=cell.group.list[[celltype_2]],
    min.pct.cells.genes=5,
    min.pct.cells.sj=5,
    min.gene.norm=1.0,
    seed=1,
    n.iterations=100,
    downsample=TRUE,
    show.progress=FALSE
)

# Compare genes
marvel <- CompareValues.Genes.10x(
    MarvelObject=marvel,
    show.progress=FALSE
)

head(marvel$DE$SJ$Table)

# Volcano plot
marvel <- PlotDEValues.SJ.10x(
    MarvelObject=marvel,
    pval=0.1,
    delta=1,
    min.gene.norm=1.0,
    anno=FALSE
)

p3 <- marvel$DE$SJ$VolcanoPlot$SJ$Plot
ggsave(p3, filename=paste0(output_celltype, "/fig3_Volcano.png"), dpi=1200, width=4, height=2.5)

head(marvel$DE$SJ$VolcanoPlot$SJ$Data[, c("coord.intron", "gene_short_name", "sig")])

# Isoform switch analysis
marvel <- IsoSwitch.10x(
    MarvelObject=marvel,
    pval.sj=0.05,
    delta.sj=5,
    min.gene.norm=1.0,
    pval.adj.gene=0.05,
    log2fc.gene=0.5
)

p4 <- marvel$SJ.Gene.Cor$Proportion$Plot
ggsave(p4, filename=paste0(output_celltype, "/fig4_IsoSwitch.pdf"), dpi=300)
marvel$SJ.Gene.Cor$Proportion$Table


# =============================================================================
# 8. Export DE gene lists
# =============================================================================
# Up-regulated genes
pos_up <- which(marvel$DE$SJ$VolcanoPlot$SJ$Data$sig == "up")
genes_up <- names(table(marvel$DE$SJ$VolcanoPlot$SJ$Data$gene_short_name[pos_up]))
write.csv(genes_up, file=paste0(output_celltype, "/", celltype_1, "-", celltype_2, "_up.csv"), quote=FALSE, row.names=FALSE)

# Down-regulated genes
pos_down <- which(marvel$DE$SJ$VolcanoPlot$SJ$Data$sig == "down")
genes_down <- names(table(marvel$DE$SJ$VolcanoPlot$SJ$Data$gene_short_name[pos_down]))
write.csv(genes_down, file=paste0(output_celltype, "/", celltype_1, "-", celltype_2, "_down.csv"), quote=FALSE, row.names=FALSE)

# Correlation data
write.csv(marvel$SJ.Gene.Cor$Data[, c("coord.intron", "gene_short_name", "cor.complete")], 
    file=paste0(output_celltype, "/", celltype_1, "-", celltype_2, "_Proportion.csv"), quote=FALSE)


# =============================================================================
# 9. Gene-specific SJ figures (curated gene list)
# =============================================================================
# Patch wiggleplotr::plotTranscripts for MARVEL compatibility (drop obsolete args)
library(wiggleplotr)

.orig_plotTranscripts <- wiggleplotr::plotTranscripts
unlockBinding("plotTranscripts", asNamespace("wiggleplotr"))
assignInNamespace(
  x = "plotTranscripts",
  value = function(
    exons,
    cdss = NULL,
    transcript_annotations = NULL,
    rescale_introns = TRUE,
    new_intron_length = 50,
    flanking_length = c(50, 50),
    connect_exons = TRUE,
    transcript_label = TRUE,
    region_coords = NULL,
    ...
  ) {
    dots <- list(...)
    dots$anno.colors <- NULL
    dots$anno.label.size <- NULL

    do.call(
      .orig_plotTranscripts,
      c(
        list(
          exons = exons,
          cdss = cdss,
          transcript_annotations = transcript_annotations,
          rescale_introns = rescale_introns,
          new_intron_length = new_intron_length,
          flanking_length = flanking_length,
          connect_exons = connect_exons,
          transcript_label = transcript_label,
          region_coords = region_coords
        ),
        dots
      )
    )
  },
  ns = "wiggleplotr"
)

lockBinding("plotTranscripts", asNamespace("wiggleplotr"))

# Curated genes (comments from original notebook):
#   mCAF-tCAF: PDIA3, CALM1, IQGAP1, ROCK1, NPM1, PTP4A2
#   mCAF-dCAF: PTMA, SULF1, MYH9, LOX, VIM
# KNOWN: default contrast above is mCAF vs tCAF, but the active list below
# follows the mCAF-dCAF note — switch as needed (see KNOWN_ISSUES.md).
# genes_of_interest <- c("PDIA3", "CALM1", "IQGAP1", "ROCK1", "NPM1", "PTP4A2")
genes_of_interest <- c("PTMA", "SULF1", "MYH9", "LOX", "VIM")

for (gene in genes_of_interest){
print(paste("Analyzing gene:", gene))
gene_dir <- paste0(output_celltype, "/", gene)
dir.create(gene_dir, showWarnings=FALSE, recursive=TRUE)
print('fig 4')
# Gene expression
marvel <- adhocGene.TabulateExpression.Gene.10x(
    MarvelObject=marvel,
    cell.group.list=cell.group.list,
    gene_short_name=gene,
    min.pct.cells=5,
    downsample=TRUE
)
marvel$adhocGene$Expression$Gene$Table
p10 <- marvel$adhocGene$Expression$Gene$Plot
ggsave(p10, filename=paste0(gene_dir, "/fig4_exp_", gene, ".pdf"), dpi=600)
print('fig 5')
# SJ expression (PSI)
marvel <- adhocGene.TabulateExpression.PSI.10x(
    MarvelObject=marvel,
    min.pct.cells=5
)
head(marvel$adhocGene$Expression$PSI$Table)
p6 <- marvel$adhocGene$Expression$PSI$Plot
ggsave(p6, filename=paste0(gene_dir, "/fig5_sj_", gene, ".pdf"), dpi=600)

# Differential analysis
marvel <- adhocGene.DE.Gene.10x(MarvelObject=marvel)
marvel <- adhocGene.DE.PSI.10x(MarvelObject=marvel)
results <- marvel$adhocGene$DE$PSI$Data
head(results)
print('fig 6')
    # Plot each SJ
    len <- length(table(results$figure.column))
    for (i in 1:len){
        print(paste0("Processing SJ-", i))
        coord.intron <- results[which(results$figure.column == paste0("SJ-", i)), "coord.intron"]
        coord.intron <- unique(coord.intron)
        
        # Volcano plot
        marvel <- adhocGene.PlotDEValues.10x(
            MarvelObject=marvel,
            coord.intron=coord.intron,
            log2fc.gene=0.5,
            delta.sj=5,
            label.size=2,
            point.size=3.5,
            xmin=-2.0,
            xmax=2.0,
            ymin=-25,
            ymax=25
        )
        p7 <- marvel$adhocGene$DE$VolcanoPlot$Plot
        ggsave(p7, filename=paste0(gene_dir, "/fig6_sj", i, "_VolcanoPlot.pdf"), dpi=600)
        print('fig 7')
        # SJ position visualization
        marvel <- adhocGene.PlotSJPosition.10x(
            MarvelObject=marvel,
            coord.intron=coord.intron,
            rescale_introns=FALSE,
            show.protein.coding.only=TRUE,
            anno.label.size=1.5
        )
        p8 <- marvel$adhocGene$SJPosition$Plot
        ggsave(p8, filename=paste0(gene_dir, "/fig7_sj", i, "_region.pdf"), dpi=600)
    }
}

# Save MARVEL object (core pipeline checkpoint)
save(marvel, file = paste0(output_base, "/marvel_CAF_", sample_suffix, ".rds"))


# =============================================================================
# 10. EXPLORATORY — auto-select candidate genes (High ΔPSI / isoform switch)
#     Note: genes_of_interest is later hard-set to "TPM2" for a focused plot.
# =============================================================================

library(dplyr)

CUTOFF_DELTA_PSI <- 10  # |Delta PSI| > 10 (%)

# --- High Delta PSI events ---
high_delta_candidates <- marvel$DE$SJ$VolcanoPlot$SJ$Data %>%
  filter(sig %in% c("up", "down")) %>%
  filter(abs(delta) > CUTOFF_DELTA_PSI) %>%
  arrange(desc(abs(delta)))

print(paste("Found", nrow(high_delta_candidates), "High Delta PSI events."))
write.csv(high_delta_candidates,
          file = paste0(output_celltype, "/Candidates_HighDeltaPSI.csv"),
          row.names = FALSE)

# --- Isoform switching candidates ---
iso_switch_candidates <- marvel$SJ.Gene.Cor$Data %>%
  filter(cor.complete != "n.s.") %>%
  filter(abs(delta) > CUTOFF_DELTA_PSI) %>%
  arrange(desc(abs(delta)))

print(paste("Found", nrow(iso_switch_candidates), "Isoform Switch candidates."))
write.csv(iso_switch_candidates,
          file = paste0(output_celltype, "/Candidates_IsoSwitch.csv"),
          row.names = FALSE)

final_targets <- unique(c(
  head(high_delta_candidates$gene_short_name, 30),
  head(iso_switch_candidates$gene_short_name, 30)
))

print("Top Priority Genes for Validation:")
print(final_targets)

# Gene-specific plots for auto-selected targets (overridden to TPM2 below)
genes_of_interest <- final_targets
genes_of_interest <- "TPM2"

if(length(genes_of_interest) == 0) {
  warning("No genes passed the threshold; consider lowering CUTOFF_DELTA_PSI (for example, to 5).")
} else {
  # Plot each retained gene.
  for (gene in genes_of_interest){
    print(paste("Analyzing gene:", gene))
    gene_dir <- paste0(output_celltype, "/", gene)
    dir.create(gene_dir, showWarnings=FALSE, recursive=TRUE)
    
    # 1. Gene Expression Plot
    print(paste("Plotting Expression for", gene))
    # Low-expression genes can fail at this step, so guard the plot with tryCatch.
    try({
        marvel <- adhocGene.TabulateExpression.Gene.10x(
            MarvelObject=marvel,
            cell.group.list=cell.group.list,
            gene_short_name=gene,
            min.pct.cells=5,
            downsample=TRUE
        )
        p10 <- marvel$adhocGene$Expression$Gene$Plot
        ggsave(p10, filename=paste0(gene_dir, "/fig4_exp_", gene, ".pdf"), dpi=600, width=5, height=4)
    })

    # 2. PSI Plot
    print(paste("Plotting PSI for", gene))
    try({
        marvel <- adhocGene.TabulateExpression.PSI.10x(
            MarvelObject=marvel,
            min.pct.cells=5
        )
        p6 <- marvel$adhocGene$Expression$PSI$Plot
        ggsave(p6, filename=paste0(gene_dir, "/fig5_sj_", gene, ".pdf"), dpi=600, width=5, height=4)
    })

    # 3. Differential Analysis & Volcano (Gene Level)
    try({
        marvel <- adhocGene.DE.Gene.10x(MarvelObject=marvel)
        marvel <- adhocGene.DE.PSI.10x(MarvelObject=marvel)
        results <- marvel$adhocGene$DE$PSI$Data
        
        # Plot each SJ
        len <- length(table(results$figure.column))
        if(len > 0) {
            for (i in 1:len){
                print(paste0("Processing SJ-", i))
                sj_id <- paste0("SJ-", i)
                # Retain only the row for this splice junction.
                current_sj_data <- results[which(results$figure.column == sj_id), ]
                coord.intron <- unique(current_sj_data$coord.intron)
                
                if(length(coord.intron) > 0) {
                    # Volcano plot specific to this gene's SJs
                    marvel <- adhocGene.PlotDEValues.10x(
                        MarvelObject=marvel,
                        coord.intron=coord.intron,
                        log2fc.gene=0.5,
                        delta.sj=5,
                        label.size=2,
                        point.size=3.5,
                        xmin=-2.0,
                        xmax=2.0,
                        ymin=-25,
                        ymax=25
                    )
                    p7 <- marvel$adhocGene$DE$VolcanoPlot$Plot
                    ggsave(p7, filename=paste0(gene_dir, "/fig6_sj", i, "_VolcanoPlot.pdf"), dpi=600)
                    
                    # SJ position visualization
                    marvel <- adhocGene.PlotSJPosition.10x(
                        MarvelObject=marvel,
                        coord.intron=coord.intron,
                        rescale_introns=FALSE,
                        show.protein.coding.only=TRUE,
                        anno.label.size=1.5
                    )
                    p8 <- marvel$adhocGene$SJPosition$Plot
                    ggsave(p8, filename=paste0(gene_dir, "/fig7_sj", i, "_region.pdf"), dpi=600, width=8, height=4)
                }
            }
        }
    })
  }
}











# =============================================================================
# 11. EXPLORATORY — require >= MIN_SIG_SJ_COUNT significant SJs per gene
#     (prefer genes with multiple robust SJ changes; plots under de_filter/)
# =============================================================================
library(dplyr)

CUTOFF_DELTA_PSI <- 10
MIN_SIG_SJ_COUNT <- 3

print(paste(
  "Filter: |Delta PSI| >", CUTOFF_DELTA_PSI,
  "and >=", MIN_SIG_SJ_COUNT, "significant SJs per gene"
))

high_delta_candidates <- marvel$DE$SJ$VolcanoPlot$SJ$Data %>%
  filter(sig %in% c("up", "down")) %>%
  filter(abs(delta) > CUTOFF_DELTA_PSI) %>%
  add_count(gene_short_name, name = "sj_count") %>%
  filter(sj_count >= MIN_SIG_SJ_COUNT) %>%
  arrange(desc(sj_count), desc(abs(delta)))

print(paste(
  "Found", length(unique(high_delta_candidates$gene_short_name)),
  "genes with >=", MIN_SIG_SJ_COUNT, "High Delta PSI SJs."
))

write.csv(
  high_delta_candidates,
  file = paste0(output_celltype, "/Candidates_HighDeltaPSI_Filtered.csv"),
  row.names = FALSE
)

iso_switch_candidates <- marvel$SJ.Gene.Cor$Data %>%
  filter(cor.complete != "n.s.") %>%
  filter(abs(delta) > CUTOFF_DELTA_PSI) %>%
  add_count(gene_short_name, name = "sj_count") %>%
  filter(sj_count >= MIN_SIG_SJ_COUNT) %>%
  arrange(desc(sj_count), desc(abs(delta)))

print(paste(
  "Found", length(unique(iso_switch_candidates$gene_short_name)),
  "genes with >=", MIN_SIG_SJ_COUNT, "Isoform Switch SJs."
))

write.csv(
  iso_switch_candidates,
  file = paste0(output_celltype, "/Candidates_IsoSwitch_Filtered.csv"),
  row.names = FALSE
)

final_targets <- unique(c(
  high_delta_candidates$gene_short_name,
  iso_switch_candidates$gene_short_name
))

print("Top Priority Genes for Validation (High SJ Count):")
print(final_targets)

genes_of_interest <- final_targets

if (length(genes_of_interest) == 0) {
  warning("No genes passed MIN_SIG_SJ_COUNT / Delta PSI filters.")
} else {
  for (gene in genes_of_interest) {
    print(paste("Analyzing gene:", gene))
    gene_dir <- paste0(output_celltype, "/de_filter/", gene)
    dir.create(gene_dir, showWarnings=FALSE, recursive=TRUE)
    
    # 1. Gene Expression Plot
    try({
        marvel <- adhocGene.TabulateExpression.Gene.10x(
            MarvelObject=marvel,
            cell.group.list=cell.group.list,
            gene_short_name=gene,
            min.pct.cells=5,
            downsample=TRUE
        )
        p10 <- marvel$adhocGene$Expression$Gene$Plot
        ggsave(p10, filename=paste0(gene_dir, "/fig4_exp_", gene, ".pdf"), dpi=600, width=5, height=4)
    })

    # 2. PSI Plot
    try({
        marvel <- adhocGene.TabulateExpression.PSI.10x(
            MarvelObject=marvel,
            min.pct.cells=5
        )
        p6 <- marvel$adhocGene$Expression$PSI$Plot
        ggsave(p6, filename=paste0(gene_dir, "/fig5_sj_", gene, ".pdf"), dpi=600, width=5, height=4)
    })

    # 3. Differential Analysis & Volcano (Gene Level)
    try({
        marvel <- adhocGene.DE.Gene.10x(MarvelObject=marvel)
        marvel <- adhocGene.DE.PSI.10x(MarvelObject=marvel)
        results <- marvel$adhocGene$DE$PSI$Data
        
        # Plot each SJ
        len <- length(table(results$figure.column))
        if(len > 0) {
            for (i in 1:len){
                # Progress output only.
                # print(paste0("Processing SJ-", i))
                
                sj_id <- paste0("SJ-", i)
                current_sj_data <- results[which(results$figure.column == sj_id), ]
                coord.intron <- unique(current_sj_data$coord.intron)
                
                if(length(coord.intron) > 0) {
                    # Volcano plot
                    marvel <- adhocGene.PlotDEValues.10x(
                        MarvelObject=marvel,
                        coord.intron=coord.intron,
                        log2fc.gene=0.5,
                        delta.sj=5,
                        label.size=2,
                        point.size=3.5,
                        xmin=-2.0,
                        xmax=2.0,
                        ymin=-25,
                        ymax=25
                    )
                    p7 <- marvel$adhocGene$DE$VolcanoPlot$Plot
                    ggsave(p7, filename=paste0(gene_dir, "/fig6_sj", i, "_VolcanoPlot.pdf"), dpi=600)
                    
                    # SJ position
                    marvel <- adhocGene.PlotSJPosition.10x(
                        MarvelObject=marvel,
                        coord.intron=coord.intron,
                        rescale_introns=FALSE,
                        show.protein.coding.only=TRUE,
                        anno.label.size=1.5
                    )
                    p8 <- marvel$adhocGene$SJPosition$Plot
                    ggsave(p8, filename=paste0(gene_dir, "/fig7_sj", i, "_region.pdf"), dpi=600, width=8, height=4)
                }
            }
        }
    })
  }
}



# =============================================================================
# 12. EXPLORATORY — short genes with multiple significant SJs (plots under len_filter/)
#     Comment thresholds disagree with variable names in the original notebook
#     (e.g. "30k" vs CUTOFF_LENGTH=50000); see KNOWN_ISSUES.md.
# =============================================================================
library(dplyr)
library(stringr)

print("Calculating gene lengths from GTF...")

get_gene_lengths <- function(marvel_obj) {
  if (is.null(marvel_obj$gtf)) stop("MARVEL object does not contain GTF data.")
  df_gtf <- marvel_obj$gtf %>%
    filter(V3 == "gene") %>%
    mutate(gene_length = V5 - V4)
  df_gtf$gene_short_name <- str_extract(df_gtf$V9, 'gene_name "([^"]+)"')
  df_gtf$gene_short_name <- gsub('gene_name "|"', "", df_gtf$gene_short_name)
  return(df_gtf[, c("gene_short_name", "gene_length")])
}

df_gene_lengths <- get_gene_lengths(marvel)
print(paste("Calculated lengths for", nrow(df_gene_lengths), "genes."))

CUTOFF_DELTA_PSI <- 3
CUTOFF_SJ_COUNT <- 2
CUTOFF_LENGTH <- 50000

df_sj_stats <- marvel$DE$SJ$VolcanoPlot$SJ$Data
delta_col <- intersect(colnames(df_sj_stats), c("delta", "mean.diff", "psi.diff"))[1]

sig_sj_counts <- df_sj_stats %>%
  filter(sig %in% c("up", "down")) %>%
  filter(abs(get(delta_col)) > CUTOFF_DELTA_PSI) %>%
  group_by(gene_short_name) %>%
  summarise(n_sig_sj = n()) %>%
  filter(n_sig_sj > CUTOFF_SJ_COUNT)

targets_short_complex <- sig_sj_counts %>%
  inner_join(df_gene_lengths, by = "gene_short_name") %>%
  filter(gene_length < CUTOFF_LENGTH) %>%
  arrange(desc(n_sig_sj))

print(paste("Found", nrow(targets_short_complex), "target genes."))
print(targets_short_complex)

write.csv(
  targets_short_complex,
  file = paste0(output_celltype, "/Targets_ShortGene_ComplexSJ.csv"),
  row.names = FALSE
)

genes_of_interest <- targets_short_complex$gene_short_name

if (length(genes_of_interest) == 0) {
  warning("No short/complex SJ genes passed filters.")
} else {
  for (gene in genes_of_interest) {
    print(paste("Analyzing gene:", gene))
    gene_dir <- paste0(output_celltype, "/len_filter/", gene)
    dir.create(gene_dir, showWarnings=FALSE, recursive=TRUE)
    
    # 1. Gene Expression Plot
    try({
        marvel <- adhocGene.TabulateExpression.Gene.10x(
            MarvelObject=marvel,
            cell.group.list=cell.group.list,
            gene_short_name=gene,
            min.pct.cells=5,
            downsample=TRUE
        )
        p10 <- marvel$adhocGene$Expression$Gene$Plot
        ggsave(p10, filename=paste0(gene_dir, "/fig4_exp_", gene, ".pdf"), dpi=600, width=5, height=4)
    })

    # 2. PSI Plot
    try({
        marvel <- adhocGene.TabulateExpression.PSI.10x(
            MarvelObject=marvel,
            min.pct.cells=5
        )
        p6 <- marvel$adhocGene$Expression$PSI$Plot
        ggsave(p6, filename=paste0(gene_dir, "/fig5_sj_", gene, ".pdf"), dpi=600, width=5, height=4)
    })

    # 3. Differential Analysis & Volcano (Gene Level)
    try({
        marvel <- adhocGene.DE.Gene.10x(MarvelObject=marvel)
        marvel <- adhocGene.DE.PSI.10x(MarvelObject=marvel)
        results <- marvel$adhocGene$DE$PSI$Data
        
        # Plot each SJ
        len <- length(table(results$figure.column))
        if(len > 0) {
            for (i in 1:len){
                # Progress output only.
                # print(paste0("Processing SJ-", i))
                
                sj_id <- paste0("SJ-", i)
                current_sj_data <- results[which(results$figure.column == sj_id), ]
                coord.intron <- unique(current_sj_data$coord.intron)
                
                if(length(coord.intron) > 0) {
                    # Volcano plot
                    marvel <- adhocGene.PlotDEValues.10x(
                        MarvelObject=marvel,
                        coord.intron=coord.intron,
                        log2fc.gene=0.5,
                        delta.sj=5,
                        label.size=2,
                        point.size=3.5,
                        xmin=-2.0,
                        xmax=2.0,
                        ymin=-25,
                        ymax=25
                    )
                    p7 <- marvel$adhocGene$DE$VolcanoPlot$Plot
                    ggsave(p7, filename=paste0(gene_dir, "/fig6_sj", i, "_VolcanoPlot.pdf"), dpi=600)
                    
                    # SJ position
                    marvel <- adhocGene.PlotSJPosition.10x(
                        MarvelObject=marvel,
                        coord.intron=coord.intron,
                        rescale_introns=FALSE,
                        show.protein.coding.only=TRUE,
                        anno.label.size=1.5
                    )
                    p8 <- marvel$adhocGene$SJPosition$Plot
                    ggsave(p8, filename=paste0(gene_dir, "/fig7_sj", i, "_region.pdf"), dpi=600, width=8, height=4)
                }
            }
        }
    })
  }
}



