#!/usr/bin/env python3
# =============================================================================
# Module : 2-CAF / Post-clustering CAF analysis + MARVEL metadata export
# Dataset: High-confidence CAF subtypes (run_paper.h5ad from 5-caf-cluster.py)
# Source : 2_hBreast/bin/old/CAF_anaylsis.py
#
# Workflow:
#   1) DEG (Wilcoxon) per CAF subtype → GeneList_celltype.csv
#   2) KEGG enrichment (gseapy enrichr) for top genes per subtype
#   3) Export meta.csv for MARVEL / SJ analysis (step 4)
#        columns: cell_id, sample, barcode, celltype
#
# Prerequisite:
#   .../CAFs/unsupervised/filter_0.7/run_paper.h5ad
#   .../BreastCancer_major/cls_total.h5ad
# =============================================================================

import os

import gseapy as gp
import matplotlib.pyplot as plt
import pandas as pd
import scanpy as sc
from gseapy.plot import barplot

plt.rcParams["pdf.fonttype"] = 42

FILTER_DIR = (
    "/path/to/floraseq/"
    "caf/result/BreastCancer_major/"
    "CAFs/unsupervised/filter_0.7"
)
CLS_TOTAL = (
    "/path/to/floraseq/"
    "caf/result/BreastCancer_major/cls_total.h5ad"
)
MARVEL_DIR = (
    "/path/to/floraseq/marvel/CAF_sc"
)


# =============================================================================
# 1. Load annotated high-confidence CAFs
# =============================================================================
os.chdir(FILTER_DIR)
adata_total = sc.read_h5ad("run_paper.h5ad")
cls = sc.read_h5ad(CLS_TOTAL)
cls = cls[adata_total.obs.index, :].copy()
cls.obs["celltype"] = adata_total.obs["celltype"]

if "log1p" in adata_total.uns and isinstance(adata_total.uns["log1p"], dict):
    adata_total.uns["log1p"]["base"] = None


# =============================================================================
# 2. Differential genes per CAF subtype
# =============================================================================
sc.tl.rank_genes_groups(adata_total, "celltype", method="wilcoxon")

celltypes = list(adata_total.obs.celltype.value_counts().index)
gs = pd.DataFrame()
for ct in celltypes:
    gs[ct] = list(adata_total.uns["rank_genes_groups"]["names"][ct][0:100])

gs.to_csv("GeneList_celltype.csv")


# =============================================================================
# 3. Enrichment analysis (top 50 genes per subtype)
# =============================================================================
# by metascope


# =============================================================================
# 4. Export barcode metadata for MARVEL / SJ (step 4)
# =============================================================================
bar = []
sample = []
merge = []
celltype = []
for i in range(cls.obs.shape[0]):
    bar.append(cls.obs.index[i].split("-")[0])
    sample.append(cls.obs["sample"][i])
    merge.append("human_" + sample[i] + "_3r_" + bar[i])
    celltype.append(cls.obs.celltype.iloc[i])

out = pd.DataFrame(merge)
out.columns = ["cell_id"]
out["sample"] = sample
out["barcode"] = bar
out.index = out["cell_id"]
out["celltype"] = celltype

os.makedirs(MARVEL_DIR, exist_ok=True)
os.chdir(MARVEL_DIR)
out.to_csv("meta.csv", sep="\t")
print(f"Saved: {os.path.join(MARVEL_DIR, 'meta.csv')}")
