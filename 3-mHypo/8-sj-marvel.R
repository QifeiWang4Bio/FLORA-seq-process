#!/usr/bin/env Rscript
# =============================================================================
# Module : 3-mHypo / MARVEL SJ on C25 neurons (GABA-1 vs GLU-2)
# Source : 3_mHypo/bin/old/MARVEL_DropBased.sh (R section)
# Prerequisite: 7-sj-process.sh; 4-anno-coarse.py → meta.csv
# Environment: conda activate sc-marvel / marvel_; Seurat, MARVEL, tidyverse
# =============================================================================

library(tidyverse)
library(Seurat)
library(R.utils)
library(Matrix)
library(MARVEL)
library(ggplot2)
library(patchwork)

setwd("/path/to/floraseq/marvel/proj_mHypo")

# FINAL: C25 neuron meta (written by 4-anno-coarse.py)
meta <- read.table(
  "/path/to/floraseq/mhypo/seurat_out/mHypo_meng_hvg5k_C25/meta.csv",
  sep = "\t", header = TRUE, row.names = 1
)

# EXPLORATORY (C66 meta) — do not use for the final track:
# meta <- read.table(
#   ".../mHypo_meng_hvg8k_filter_C66_bin/meta.csv",
#   sep = "\t", header = TRUE, row.names = 1
# )

samples <- paste0("group", seq(1, 24), "_3r")

# =============================================================================
# 1) Raw gene count (subset to meta barcodes)
# =============================================================================
dir.create("./gene", showWarnings = FALSE)
dir.create("./gene/raw", showWarnings = FALSE)
dir.create("./gene/norm", showWarnings = FALSE)

sce.list <- lapply(samples, function(x) {
  cat(paste0("\n", x))
  fls <- list.files(paste0("./star/", x, "/Solo.out/Gene/filtered"), full.name = TRUE)
  if (length(grep("gz", fls)) != 3) {
    for (fl in fls) { gzip(fl, remove = FALSE) }
  }
  count <- Read10X(paste0("./star/", x, "/Solo.out/Gene/filtered"))
  sce <- CreateSeuratObject(count)
  sce <- RenameCells(sce, add.cell.id = x)
  sce$Group <- x
  bar <- intersect(colnames(sce), rownames(meta))
  sce <- sce[, bar]
  print(paste("Cell counts:", dim(sce)[2]))
  return(sce)
})
names(sce.list) <- samples

sce <- merge(sce.list[[1]], sce.list[-1])
table(sce$Group)

count_merge <- sce@assays$RNA@counts
writeMM(count_merge, file = "./gene/raw/matrix.mtx")

count_merge_phe <- sce@meta.data %>%
  tibble::rownames_to_column("cell.id") %>%
  dplyr::select(cell.id)
write.table(count_merge_phe, sep = "\t", row.names = FALSE,
            file = "./gene/raw/phenoData.txt")

count_merge_feat <- data.frame(gene_short_name = rownames(sce))
write.table(count_merge_feat, sep = "\t", row.names = FALSE,
            file = "./gene/raw/featureData.txt")

# =============================================================================
# 2) Normalised gene count
# =============================================================================
sce <- NormalizeData(sce)
count_merge_norm <- Matrix::Matrix(exp(sce@assays$RNA@data) - 1, sparse = TRUE)
writeMM(count_merge_norm, file = "./gene/norm/matrix_normalised.mtx")

count_merge_norm_phe <- sce@meta.data %>%
  tibble::rownames_to_column("cell.id") %>%
  dplyr::select(cell.id, Group)
write.table(count_merge_norm_phe, sep = "\t", row.names = FALSE,
            file = "./gene/norm/phenoData.txt")

count_merge_norm_feat <- data.frame(gene_short_name = rownames(sce))
write.table(count_merge_norm_feat, sep = "\t", row.names = FALSE,
            file = "./gene/norm/featureData.txt")

# =============================================================================
# 3) Raw SJ count
# =============================================================================
dir.create("./sj/raw", showWarnings = FALSE, recursive = TRUE)

for (x in samples) {
  print(x)
  sj_name <- data.table::fread(paste0("./star/", x, "/Solo.out/SJ.out.tab"), data.table = FALSE) %>%
    dplyr::mutate(coord.intron = paste(V1, V2, V3, sep = ":")) %>%
    dplyr::select(coord.intron)
  write.table(sj_name, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE,
              file = paste0("./star/", x, "/Solo.out/SJ/raw/features.tsv"))
}

sce_sj.list <- lapply(samples, function(x) {
  print(x)
  fls <- list.files(paste0("./star/", x, "/Solo.out/SJ/raw"), full.name = TRUE)
  if (length(grep("gz", fls)) != 3) {
    for (fl in fls) { gzip(fl, remove = FALSE) }
  }
  count <- Read10X(paste0("./star/", x, "/Solo.out/SJ/raw"), gene.column = 1)
  sce_sj <- CreateSeuratObject(count)
  sce_sj <- RenameCells(sce_sj, add.cell.id = x)
  sce_sj <- sce_sj[, colnames(sce.list[[x]])]
  return(sce_sj)
})

sce_sj <- merge(sce_sj.list[[1]], sce_sj.list[-1])
sj_count_merge <- sce_sj@assays$RNA@counts
writeMM(sj_count_merge, file = "./sj/raw/matrix.mtx")

sj_count_merge_phe <- sce_sj@meta.data %>%
  tibble::rownames_to_column("cell.id") %>%
  dplyr::select(cell.id)
write.table(sj_count_merge_phe, sep = "\t", row.names = FALSE,
            file = "./sj/raw/phenoData.txt")

sj_count_merge_feat <- data.frame(coord.intron = rownames(sce_sj)) %>%
  dplyr::mutate(coord.intron = gsub("-", "_", coord.intron))
write.table(sj_count_merge_feat, sep = "\t", row.names = FALSE,
            file = "./sj/raw/featureData.txt")

# =============================================================================
# 4) UMAP coordinates for MARVEL
# =============================================================================
dir.create("./umap", showWarnings = FALSE)

sce <- sce %>%
  NormalizeData() %>%
  FindVariableFeatures() %>%
  ScaleData() %>%
  Seurat::RunPCA() %>%
  Seurat::RunUMAP(dims = 1:30)

umap_embed <- sce@reductions$umap@cell.embeddings %>%
  as.data.frame() %>%
  tibble::rownames_to_column("cell.id") %>%
  dplyr::rename(x = UMAP_1, y = UMAP_2)
write.table(umap_embed, row.names = FALSE, sep = "\t",
            file = "./umap/UMAP_coordinate.txt")

# =============================================================================
# 5) Create MARVEL object
# =============================================================================
df.gene.norm <- readMM("./gene/norm/matrix_normalised.mtx")
df.gene.norm.pheno <- read.table("./gene/norm/phenoData.txt", sep = "\t", header = TRUE)
df.gene.norm.feature <- read.table("./gene/norm/featureData.txt", sep = "\t", header = TRUE)

df.gene.count <- readMM("./gene/raw/matrix.mtx")
df.gene.count.pheno <- read.table("./gene/raw/phenoData.txt", sep = "\t", header = TRUE)
df.gene.count.feature <- read.table("./gene/raw/featureData.txt", sep = "\t", header = TRUE)

df.sj.count <- readMM("./sj/raw/matrix.mtx")
df.sj.count.pheno <- read.table("./sj/raw/phenoData.txt", sep = "\t", header = TRUE)
df.sj.count.feature <- read.table("./sj/raw/featureData.txt", sep = "\t", header = TRUE)
row.names(df.sj.count) <- df.sj.count.feature$coord.intron
row.names(df.sj.count) <- paste("chr", row.names(df.sj.count), sep = "")
df.sj.count.feature$coord.intron <- row.names(df.sj.count)

rownames(df.sj.count) <- df.sj.count.feature[, 1]
colnames(df.sj.count) <- df.sj.count.pheno[, 1]

df.coord <- read.table("./umap/UMAP_coordinate.txt", sep = "\t", header = TRUE)

# Align meta to SJ columns (C25)
meta <- meta[colnames(df.sj.count), ]

gtf <- data.table::fread("../basic/gtf/genes.gtf", sep = "\t", header = FALSE) %>%
  as.data.frame()
gtf$V9 <- gsub("gene_type", "gene_biotype", gtf$V9)
gtf$V1 <- gsub("chr", "", gtf$V1)

marvel <- CreateMarvelObject.10x(
  gene.norm.matrix = df.gene.norm,
  gene.norm.pheno = df.gene.norm.pheno,
  gene.norm.feature = df.gene.norm.feature,
  gene.count.matrix = df.gene.count,
  gene.count.pheno = df.gene.count.pheno,
  gene.count.feature = df.gene.count.feature,
  sj.count.matrix = df.sj.count,
  sj.count.pheno = df.sj.count.pheno,
  sj.count.feature = df.sj.count.feature,
  pca = df.coord,
  gtf = gtf
)

# =============================================================================
# 6) Annotate / filter SJ
# =============================================================================
marvel <- AnnotateGenes.10x(MarvelObject = marvel)
marvel <- AnnotateSJ.10x(MarvelObject = marvel)
marvel <- ValidateSJ.10x(MarvelObject = marvel)
marvel <- FilterGenes.10x(MarvelObject = marvel, gene.type = "protein_coding")

rownames(marvel$sample.metadata) <- marvel$sample.metadata$cell.id
# FINAL: C25 labels as written in meta.csv (e.g. GABA-1, GLU-2)
marvel$sample.metadata$celltype <- meta$celltype

# EXPLORATORY (C66 label cleanup) — only needed for C66 meta:
# meta$celltype <- gsub("-", "_", meta$celltype)
# marvel$sample.metadata$celltype <- meta$celltype

cell.group.list <- split(marvel$sample.metadata$cell.id, marvel$sample.metadata$celltype)

# =============================================================================
# 7) FINAL contrast: GABA-1 vs GLU-2 (C25)
# =============================================================================
dir.create("C25", showWarnings = FALSE)
celltype_1 <- "GABA-1"
celltype_2 <- "GLU-2"

marvel <- PlotPctExprCells.Genes.10x(
  MarvelObject = marvel,
  cell.group.g1 = cell.group.list[[celltype_1]],
  cell.group.g2 = cell.group.list[[celltype_2]],
  min.pct.cells = 5
)
ggsave(marvel$pct.cells.expr$Gene$Plot, filename = "C25/fig1_exp.png")

marvel <- PlotPctExprCells.SJ.10x(
  MarvelObject = marvel,
  cell.group.g1 = cell.group.list[[celltype_1]],
  cell.group.g2 = cell.group.list[[celltype_2]],
  min.pct.cells.genes = 5,
  min.pct.cells.sj = 5,
  downsample = TRUE,
  downsample.pct.sj = 100
)
ggsave(marvel$pct.cells.expr$SJ$Plot, filename = "C25/fig2_sj.png")

marvel <- CompareValues.SJ.10x(
  MarvelObject = marvel,
  cell.group.g1 = cell.group.list[[celltype_1]],
  cell.group.g2 = cell.group.list[[celltype_2]],
  min.pct.cells.genes = 5,
  min.pct.cells.sj = 5,
  min.gene.norm = 1.0,
  seed = 1,
  n.iterations = 100,
  downsample = TRUE,
  show.progress = FALSE
)
marvel <- CompareValues.Genes.10x(MarvelObject = marvel, show.progress = FALSE)

marvel <- PlotDEValues.SJ.10x(
  MarvelObject = marvel,
  pval = 0.05,
  delta = 1,
  min.gene.norm = 1.0,
  anno = FALSE
)
ggsave(marvel$DE$SJ$VolcanoPlot$SJ$Plot, filename = "C25/fig3_Volcano.png")

marvel <- IsoSwitch.10x(
  MarvelObject = marvel,
  pval.sj = 0.05,
  delta.sj = 5,
  min.gene.norm = 1.0,
  pval.adj.gene = 0.05,
  log2fc.gene = 0.5
)

# =============================================================================
# 8) Gene-specific plots (C25 candidate genes)
# =============================================================================
# C25: GABA-1_GLU-2 — Srsf11, Luc7l3, Kmt2a, Sltm, Arglu1
genes_of_interest <- c("Srsf11", "Luc7l3", "Kmt2a", "Sltm", "Arglu1")

for (gene in genes_of_interest) {
  message("Analyzing gene: ", gene)
  marvel <- adhocGene.TabulateExpression.Gene.10x(
    MarvelObject = marvel,
    cell.group.list = cell.group.list,
    gene_short_name = gene,
    min.pct.cells = 5,
    downsample = TRUE
  )
  ggsave(marvel$adhocGene$Expression$Gene$Plot,
         filename = paste0("C25/fig4_exp_", gene, ".pdf"), dpi = 600)

  marvel <- adhocGene.TabulateExpression.PSI.10x(
    MarvelObject = marvel,
    min.pct.cells = 5
  )
  ggsave(marvel$adhocGene$Expression$PSI$Plot,
         filename = paste0("C25/fig5_sj_", gene, ".pdf"), dpi = 600)
}

save(marvel, file = "marvel_C25.rds")

# =============================================================================
# EXPLORATORY — C66 subtype contrasts (from original DropBased notebook)
# Not part of the final C25 SJ track. Uncomment + switch meta to C66 to use.
# =============================================================================
#: '
# pos <- which(marvel$sample.metadata$celltype %in%
#   c("Nkx2_4.GABA_3", "Tcf4.GLU_3", "Pmch.GLU_7", "Foxb1.GLU_2"))
# cell.group.list <- split(
#   marvel$sample.metadata$cell.id[pos],
#   marvel$sample.metadata$celltype[pos]
# )
# # Example DE: Nkx2_4.GABA_3 vs Tcf4.GLU_3 → figures under C66/
# marvel <- CompareValues.SJ.10x(
#   MarvelObject = marvel,
#   cell.group.g1 = cell.group.list$"Nkx2_4.GABA_3",
#   cell.group.g2 = cell.group.list$"Tcf4.GLU_3",
#   min.pct.cells.genes = 5, min.pct.cells.sj = 5,
#   min.gene.norm = 1.0, seed = 1, n.iterations = 100,
#   downsample = TRUE, show.progress = FALSE
# )
# '

# Optional bamCoverage track (original notebook footnote):
# samtools index -@ 12 GLU-2.Aligned.sortedByCoord.out.bam
# bamCoverage -b GLU-2.Aligned.sortedByCoord.out.bam -o GLU-2.bw \
#   --numberOfProcessors 12 --binSize 10 --normalizeUsing RPGC \
#   --effectiveGenomeSize 2652783500
