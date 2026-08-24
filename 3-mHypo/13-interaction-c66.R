#!/usr/bin/env Rscript
# =============================================================================
# Module : 3-mHypo / NeuronChat cell–cell interaction at C66
# Source : 3_mHypo/bin/old/Interaction_C66_new.R
#          (+ MTX→RDS from 6-fig6-interaction.R §2)
#
# Prerequisite: python 12-interaction-prep.py
# Default panel: Metabolism & Feeding → interaction/neu_group1_metabolism/
#
# Environment: conda activate sc-chat; NeuronChat, CellChat, Seurat, ggalluvial
# =============================================================================

library(Seurat)
library(Matrix)
library(dplyr)
library(NeuronChat)
library(CellChat)
library(ggplot2)
library(ggalluvial)

C66_DIR <- "/path/to/floraseq/mhypo/seurat_out/mHypo_meng_hvg8k_filter_C66_bin"
INTER_DIR <- file.path(C66_DIR, "interaction")
RDS_FILE <- file.path(INTER_DIR, "adata_neuron_norm.rds")

close_all_devices <- function() {
  while (dev.cur() > 1) {
    dev.off()
  }
  cat("All graphics devices closed.\n")
}

# =============================================================================
# 1) MTX → Seurat RDS  (skip if RDS already exists)
# =============================================================================
setwd(INTER_DIR)

if (!file.exists(RDS_FILE)) {
  cat("Reading MTX format expression matrix...\n")
  expr_matrix <- readMM(file.path(INTER_DIR, "matrix.mtx"))
  genes <- read.table(file.path(INTER_DIR, "genes.tsv"),
                      header = FALSE, stringsAsFactors = FALSE)$V1
  barcodes <- read.table(file.path(INTER_DIR, "barcodes.tsv"),
                         header = FALSE, stringsAsFactors = FALSE)$V1

  cat(sprintf("Matrix dimensions: %d x %d\n", nrow(expr_matrix), ncol(expr_matrix)))
  cat(sprintf("Genes: %d, Barcodes: %d\n", length(genes), length(barcodes)))
  if (nrow(expr_matrix) != length(genes) || ncol(expr_matrix) != length(barcodes)) {
    stop("Dimension mismatch! Matrix should be genes x cells format.")
  }

  rownames(expr_matrix) <- genes
  colnames(expr_matrix) <- barcodes

  metadata <- read.csv(file.path(INTER_DIR, "metadata.csv"),
                       row.names = 1, check.names = FALSE)
  metadata <- metadata[colnames(expr_matrix), ]

  seurat_obj <- CreateSeuratObject(
    counts = expr_matrix,
    meta.data = metadata,
    project = "NeuronChat"
  )
  # Expression from prep is already normalized (log1p); store in data layer
  seurat_obj <- SetAssayData(seurat_obj, layer = "data", new.data = expr_matrix)

  saveRDS(seurat_obj, RDS_FILE)
  cat("Wrote", RDS_FILE, "\n")
  rm(seurat_obj, expr_matrix)
  gc()
}

# =============================================================================
# 2) NeuronChat  (Interaction_C66_new.R)
# =============================================================================
data <- readRDS(RDS_FILE)

# Harmonize grouping column
if (!("anno_l3" %in% colnames(data@meta.data))) {
  if ("C66" %in% colnames(data@meta.data)) {
    data$anno_l3 <- data$C66
  } else if ("Pred" %in% colnames(data@meta.data)) {
    data$anno_l3 <- data$Pred
  } else {
    stop("Need anno_l3 / C66 / Pred in Seurat meta.data")
  }
} else if ("C66" %in% colnames(data@meta.data) && all(is.na(data$anno_l3))) {
  data$anno_l3 <- data$C66
}

cat("Available cell types in anno_l3:\n")
print(table(data$anno_l3))

# ===== OPTION: analysis mode =====
USE_ALL_NEURONS <- FALSE

# Active group label → also used as output folder name
# Groups from Interaction_C66_new.R (smaller panels). For the longer
# fig6-interaction panels see comments below.
GROUP_NAME <- "neu_group"

cate=c('C66-6: Samd3.GLU-2','C66-2: Irx5.GLU-2','C66-14: Tcf4.GLU-3',
	'C66-21: Bace2.GLU-5','C66-17: Lpar1.GLU-4')


if (USE_ALL_NEURONS) {
  cat("\n==> Using ALL neuron cell types for analysis\n")
  obj <- "neu_all"
} else {
  cat("\n==> Using selected cell types for analysis:\n")
  print(cate)
  data <- subset(x = data, subset = anno_l3 %in% cate)
  obj <- GROUP_NAME
}

data$anno_l3 <- factor(data$anno_l3)
data <- droplevels(data)

cat(sprintf("\nAnalyzing %d cells from %d cell types\n",
            ncol(data), nlevels(data$anno_l3)))

# Create NeuronChat object
expr_mat <- tryCatch(
  GetAssayData(data, assay = "RNA", layer = "data"),
  error = function(e) GetAssayData(data, assay = "RNA", slot = "data")
)
x <- createNeuronChat(as.matrix(expr_mat), DB = "mouse", group.by = data$anno_l3)

cat("Running NeuronChat analysis...\n")
x <- run_NeuronChat(x, M = 500)
net_aggregated_x <- net_aggregation(x@net, method = "weight")

cat("\nDetected interactions:\n")
print(names(which(lapply(x@net, sum) > 0)))

dir.create(obj, showWarnings = FALSE)
while (dev.cur() > 1) dev.off()

# Aggregated visuals
tryCatch({
  png(filename = paste0(obj, "/circle_neuron.png"),
      width = 8, height = 6, units = "in", res = 600)
  netVisual_circle_neuron(
    net_aggregated_x,
    arrow.width = 0.1, arrow.size = 0.2,
    vertex.label.cex = 0.5, edge.width.max = 1
  )
  dev.off()
}, error = function(e) {
  if (dev.cur() > 1) dev.off()
  cat("Warning: Could not generate aggregated circle plot\n")
})

tryCatch({
  pdf(paste0(obj, "/chord_neuron.pdf"), width = 10, height = 10)
  netVisual_chord_neuron(x, method = "weight", lab.cex = 1)
  dev.off()
}, error = function(e) {
  if (dev.cur() > 1) dev.off()
  cat("Warning: Could not generate aggregated chord plot\n")
})

tryCatch({
  pdf(paste0(obj, "/heatmap_neuron.pdf"), width = 16, height = 10)
  heatmap_aggregated(x, method = "weight")
  dev.off()
}, error = function(e) {
  if (dev.cur() > 1) dev.off()
  cat("Warning: Could not generate aggregated heatmap\n")
})

# Per-interaction plots
cat("\nGenerating plots for individual interactions...\n")
while (dev.cur() > 1) dev.off()

for (inter_obj in names(which(lapply(x@net, sum) > 0))) {
  cat(sprintf("Plotting: %s\n", inter_obj))

  tryCatch({
    png(filename = paste0(obj, "/circle_neuron_", inter_obj, ".png"),
        width = 10, height = 10, units = "in", res = 600)
    netVisual_circle_neuron(
      x@net[[inter_obj]],
      arrow.width = 0.1, arrow.size = 0.2,
      vertex.label.cex = 0.5, edge.width.max = 1
    )
    dev.off()
  }, error = function(e) {
    if (dev.cur() > 1) dev.off()
    cat(sprintf("Warning: Could not generate circle plot for %s\n", inter_obj))
  })

  tryCatch({
    pdf(paste0(obj, "/chord_neuron_", inter_obj, ".pdf"), width = 10, height = 10)
    netVisual_chord_neuron(x, method = "weight",
                           interaction_use = inter_obj, lab.cex = 1)
    dev.off()
  }, error = function(e) {
    if (dev.cur() > 1) dev.off()
    cat(sprintf("Warning: Could not generate chord plot for %s\n", inter_obj))
  })

  tryCatch({
    png(filename = paste0(obj, "/interaction_pair_", inter_obj, ".png"),
        width = 25, height = 5, units = "in", res = 300)
    lig_tar_heatmap(x, interaction_name = inter_obj,
                    width.vector = c(0.38, 0.35, 0.27))
    dev.off()
  }, error = function(e) {
    if (dev.cur() > 1) dev.off()
    cat(sprintf("Warning: Could not generate ligand-target heatmap for %s\n", inter_obj))
  })
}

# Communication patterns
cat("\nAnalyzing communication patterns...\n")
while (dev.cur() > 1) dev.off()

tryCatch({
  pdf(paste0(obj, "/Neuron_outgoing.pdf"), width = 10, height = 10)
  x <- identifyCommunicationPatterns_Neuron(
    x, slot.name = "net", pattern = c("outgoing"), k = 4
  )
  dev.off()
}, error = function(e) {
  if (dev.cur() > 1) dev.off()
  cat("Warning: Could not generate outgoing patterns\n")
})

tryCatch({
  png(filename = paste0(obj, "/Neuron_incoming.png"),
      width = 10, height = 14, units = "in", res = 600)
  x <- identifyCommunicationPatterns_Neuron(
    x, slot.name = "net", pattern = c("incoming"), k = 4, height = 18
  )
  dev.off()
}, error = function(e) {
  if (dev.cur() > 1) dev.off()
  cat("Warning: Could not generate incoming patterns\n")
})

tryCatch({
  png(filename = paste0(obj, "/netAnalysis_river_outgoing.png"),
      width = 10, height = 14, units = "in", res = 600)
  netAnalysis_river_Neuron(
    x, slot.name = "net", pattern = c("outgoing"),
    font.size = 2.5, cutoff.1 = 0.5, cutoff.2 = 0.5
  )
  dev.off()
}, error = function(e) {
  if (dev.cur() > 1) dev.off()
  cat("Warning: Could not generate outgoing river plot\n")
})

tryCatch({
  png(filename = paste0(obj, "/netAnalysis_river_incoming.png"),
      width = 10, height = 14, units = "in", res = 600)
  netAnalysis_river_Neuron(
    x, slot.name = "net", pattern = c("incoming"),
    font.size = 2.5, cutoff.1 = 0.5, cutoff.2 = 0.5
  )
  dev.off()
}, error = function(e) {
  if (dev.cur() > 1) dev.off()
  cat("Warning: Could not generate incoming river plot\n")
})

tryCatch({
  png(filename = paste0(obj, "/rankNet_Neuron.png"),
      width = 5, height = 10, units = "in", res = 300)
  g1 <- rankNet_Neuron(x, slot.name = "net", measure = c("weight"),
                       mode = "single", font.size = 5)
  g2 <- rankNet_Neuron(x, slot.name = "net", measure = c("count"),
                       mode = "single", font.size = 5)
  print(g1 + g2)
  dev.off()
}, error = function(e) {
  if (dev.cur() > 1) dev.off()
  cat("Warning: Could not generate rankNet plot\n")
})

while (dev.cur() > 1) dev.off()
cat("\nAnalysis completed! Results saved to:", file.path(INTER_DIR, obj), "\n")
