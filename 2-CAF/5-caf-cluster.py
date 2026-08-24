#!/usr/bin/env python3
# =============================================================================
# Module : 2-CAF / Unsupervised clustering of high-confidence CAF cells
# Dataset: Breast cancer FLORA-seq CAFs
# Source : 2_hBreast/bin/old/CAFs_CellMemory.py
#
# This is step 3 of the CAF analysis track.
# =============================================================================

import os

import matplotlib.pyplot as plt
import scanpy as sc
from pandas.api.types import CategoricalDtype

plt.rcParams["pdf.fonttype"] = 42

MAJOR_DIR = (
    "/path/to/floraseq/"
    "caf/result/BreastCancer_major"
)
CAF_DIR = f"{MAJOR_DIR}/CAFs"
CLS_TOTAL = f"{MAJOR_DIR}/cls_total.h5ad"
FILTER_DIR = f"{CAF_DIR}/unsupervised/filter_0.7"


# =============================================================================
# A. EXPLORATORY — unsupervised CAF clustering WITHOUT Pred_prob filter
#    (early draft; not the final paper path)
# =============================================================================
# os.chdir(CAF_DIR)
# adata = sc.read_h5ad(f"{MAJOR_DIR}/CAFs_TotalGene.h5ad")
# ... Leiden without Pred_prob filter / draft labels ...


# =============================================================================
# B. FINAL — Pred_prob > 0.7 filter + Leiden + subtype labels (2024-04-16)
# =============================================================================
os.chdir(CAF_DIR)

# Full-gene CAF matrix (normalized/log1p from 3-prepare-annotate.py)
adata = sc.read_h5ad(f"{MAJOR_DIR}/CAFs_TotalGene.h5ad")

# Align with major-type CellMemory scores
cls = sc.read_h5ad(CLS_TOTAL)
cls = cls[adata.obs.index, :]

os.makedirs(FILTER_DIR, exist_ok=True)
os.chdir(FILTER_DIR)

adata.obs["Pred_prob"] = cls.obs["Pred_prob"]
# Optional: quick look before filtering (needs an existing embedding; skip if none)
# sc.pl.umap(adata, color="Pred_prob", frameon=False, show=False, size=20)
# plt.savefig("umap_Pred_prob_prefilter.png", dpi=300, bbox_inches="tight")

# High-confidence CAF cells only
adata = adata[adata.obs.Pred_prob > 0.7, :].copy()

if "log1p" in adata.uns and isinstance(adata.uns["log1p"], dict):
    adata.uns["log1p"]["base"] = None

# Cluster on HVG subset; keep full-gene adata for markers / DEG / export
sc.pp.highly_variable_genes(adata, n_top_genes=3000)
adata_hvg = adata[:, adata.var.highly_variable].copy()

sc.tl.pca(adata_hvg)
sc.pp.neighbors(adata_hvg, n_neighbors=20)
sc.tl.umap(adata_hvg)
sc.tl.leiden(adata_hvg)

sc.pl.umap(adata_hvg, color="leiden", frameon=False, show=False, size=20)
plt.savefig("umap_hvg3k_leiden.png", dpi=300, bbox_inches="tight")

adata.obsm["X_umap"] = adata_hvg.obsm["X_umap"]
adata.obs["leiden"] = adata_hvg.obs["leiden"]

sc.tl.rank_genes_groups(adata, "leiden", method="wilcoxon")
sc.pl.rank_genes_groups(adata, n_genes=15, sharey=False)
plt.savefig("rank_markers_hvg3k.pdf", dpi=300, bbox_inches="tight")

adata.obs["celltype"] = list(adata.obs["leiden"])
adata.obs.celltype.loc[adata.obs.celltype.isin(["0", "6", "1", "5"])] = "mCAF"
adata.obs.celltype.loc[adata.obs.celltype.isin(["3", "8"])] = "iCAF"
adata.obs.celltype.loc[adata.obs.celltype.isin(["7"])] = "tCAF"
adata.obs.celltype.loc[adata.obs.celltype.isin(["9"])] = "dCAF"
adata.obs.celltype.loc[adata.obs.celltype.isin(["2", "4"])] = "vCAF"

adata.uns["celltype_colors"] = [
    "#b15928",
    "#33a02c",
    "#a6cee3",
    "#6a3d9a",
    "#fdbf6f",
]
sc.pl.umap(adata, color="celltype", frameon=False, show=False, size=30)
plt.title("Celltype", fontsize=20)
plt.savefig("umap_hvg3k_celltype_new.png", dpi=300, bbox_inches="tight")

cat_type = CategoricalDtype(
    categories=["dCAF", "iCAF", "mCAF", "tCAF", "vCAF"], ordered=True
)
adata.obs.celltype = adata.obs.celltype.astype(cat_type)

marker_final = [
    "MKI67", "STMN1", "PTTG1",
    "C3", "FBLN1", "SELENOP",
    "COL1A1", "POSTN", "COL3A1",
    "PGK1", "TMEM158", "HSP90AA1",
    "ACTA2", "ADIRF", "SPARCL1",
]
sc.pl.matrixplot(
    adata,
    marker_final,
    "celltype",
    colorbar_title="Expression \nscaled by column",
    standard_scale="var",
    swap_axes=False,
    vmin=-1,
    vmax=1.5,
    cmap="RdBu_r",
)
plt.savefig("Heatmap_celltype.png", dpi=600, bbox_inches="tight")

for gene in ["POSTN", "HSP90AA1", "MKI67", "ACTA2", "FN1", "VIM"]:
    sc.pl.umap(adata, color=gene, frameon=False, show=False, size=30)
    plt.title(gene, fontsize=20)
    plt.savefig(f"umap_hvg3k_{gene}.png", dpi=300, bbox_inches="tight")

sc.pl.umap(adata, color="Pred_prob", frameon=False, show=False, size=20)
plt.savefig("umap_Pred_prob.png", dpi=300, bbox_inches="tight")

# Expected by 6-caf-analysis.py
adata.write("run_paper.h5ad")


