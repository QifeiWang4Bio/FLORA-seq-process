#!/usr/bin/env Rscript
# =============================================================================
# Module : 3-mHypo / QC + SoupX + DoubletFinder (24 samples)
# Source : 3_mHypo/bin/old/Seurat_rmDoublet.R
# Outputs: counts.mtx, barcode.csv, features.csv, data_filter.Rds
# Environment: conda activate myR; SoupX, Seurat, DoubletFinder
# =============================================================================

# conda activate myR
library(SoupX)
library(Seurat)
library(DoubletFinder)

# library(DropletUtils)

# 1. Load the dataset
sample <- seq(1, 24)
# sample = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 13, 17, 18, 19, 20, 21, 23)
for (i in sample){
    k=i
    print(paste0('Running:', k))
    Dir = paste0('/path/to/floraseq/mhypo/cellranger/group',k,'_3r/')
    setwd(Dir)
    temp <- load10X('outs')
    temp <- autoEstCont(temp)
    temp <- setContaminationFraction(temp, 0.2)
    temp <- adjustCounts(temp, roundToInt=TRUE)
    temp_ <- CreateSeuratObject(counts = temp, project = paste0("output_",k), min.features = 500)
    assign(paste0('output_',k), temp_)
}

data <- merge(output_1,y=c(output_2,output_3,output_4,output_5,output_6,output_7,output_8,output_9,output_10,
    output_11,output_12,output_13,output_14,output_15,output_16,output_17,output_18,output_19,output_20,
    output_21,output_22,output_23,output_24))
setwd('/path/to/floraseq/mhypo/seurat_out')

# 2. QC, mitochondrial fraction calculation and normalization.
data[["percent.mt"]] <- PercentageFeatureSet(data, pattern = "^mT-")
group <- replicate(dim(data)[2], 'mouse_hypo')
data$dataset <- group
p <- VlnPlot(data, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, group.by='dataset')
ggsave('VlnPlot_raw.png', p)

data <- subset(data, subset = nFeature_RNA > 500 & nFeature_RNA < 8000 & nCount_RNA<15000 & percent.mt < 5) # 49831 > 49232
data <- NormalizeData(data)
data <- FindVariableFeatures(data, nfeatures=3000)
data <- ScaleData(data)

# 3. find doublet
Find_doublet <- function(data){
    sweep.res.list <- paramSweep_v3(data, PCs = 1:30, sct = FALSE)
    sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
    bcmvn <- find.pK(sweep.stats)
    nExp_poi <- round(0.05*ncol(data))
    p<-as.numeric(as.vector(bcmvn[bcmvn$MeanBC==max(bcmvn$MeanBC),]$pK))
    data <- doubletFinder_v3(data, PCs = 1:30, pN = 0.25, pK = p, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
    colnames(data@meta.data)[ncol(data@meta.data)] = "doublet_info"
    return(data)
}

data <- RunPCA(data, features = VariableFeatures(object = data))
data <- Find_doublet(data)    # 25224  >  23963
data <- subset(data, subset=doublet_info=='Singlet')  # 49232 > 46770

# 4. UMAP
data <- RunPCA(data, features = VariableFeatures(object = data))
data <- FindNeighbors(data, dims = 1:30)
data <- FindClusters(data, resolution = 0.1)

data <- RunUMAP(data, dims = 1:30)
p <- DimPlot(data, reduction = "umap")
ggsave('UMAP_transfer_0.2.png', p)

saveRDS(data, 'data_filter.Rds')

# 0. Export the Seurat object as an MTX matrix.
library(Matrix)
Matrix::writeMM(t(data@assays$RNA@counts),'counts.mtx')
write.csv(rownames(data), 'features.csv', row.names=FALSE, quote=FALSE)
write.csv(colnames(data), 'barcode.csv', row.names=FALSE, quote=FALSE)


