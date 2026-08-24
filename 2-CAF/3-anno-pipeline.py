#!/usr/bin/env python3
# =============================================================================
# Module : 2-CAF / Prepare query matrix + CellMemory major annotation + CAF filter
# Dataset: Breast cancer FLORA-seq (human_sample1/2/3_3r)
# Source : 2_hBreast/bin/old/breast_CellLine.py
#
# This is step 2 of the CAF analysis track:
#   1) QC three Cell Ranger samples → write adata_caf.h5ad (post-QC counts)
#   2) Intersect genes with breast reference → BreastCancer_train/test.h5ad
#   3) CellMemory train/predict on celltype_major (legacy API block below)
#      Prefer the cleaner re-run in 4-cellmemory-total.py for cls_total.h5ad
#   4) Filter Pred == 'CAFs' from adata_caf.h5ad → CAFs_TotalGene.h5ad
#
# Final used outputs for downstream unsupervised CAF clustering (5-caf-cluster.py):
#   adata_caf.h5ad              # post-QC count matrix (all genes after MT removal)
#   BreastCancer_train.h5ad
#   BreastCancer_test.h5ad
#   CAFs_TotalGene.h5ad   (+ CAFs_test.h5ad if written separately)
#   cls_total.h5ad        (Pred, Pred_prob; typically from 4-cellmemory-total.py)
# =============================================================================

import os

import matplotlib.pyplot as plt
import numpy as np
import scanpy as sc

plt.rcParams["pdf.fonttype"] = 42

# Working / result directories from the original analysis
RESULT_DIR = (
    "/path/to/floraseq/"
    "caf/result"
)
PAPER_DIR = (
    "/path/to/floraseq/"
    "caf"
)
REF_H5AD = (
    "/path/to/external-data/cancer/breast/data.h5ad"
)


def read_data(path):
    data = sc.read_10x_h5(path)
    data.var_names_make_unique()
    data.var["mt"] = data.var_names.str.startswith("MT-")
    sc.pp.calculate_qc_metrics(
        data, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True
    )
    return data


# =============================================================================
# 1. Load three samples + QC
# =============================================================================
data_1 = read_data(
    f"{PAPER_DIR}/human_sample1_3r/outs/filtered_feature_bc_matrix.h5"
)
data_2 = read_data(
    f"{PAPER_DIR}/human_sample2_3r/outs/filtered_feature_bc_matrix.h5"
)
data_3 = read_data(
    f"{PAPER_DIR}/human_sample3_3r/outs/filtered_feature_bc_matrix.h5"
)
data_1.obs["sample"] = "sample1"
data_2.obs["sample"] = "sample2"
data_3.obs["sample"] = "sample3"

adata = sc.AnnData.concatenate(data_1, data_2, data_3)

adata.var["mt"] = adata.var_names.str.startswith("MT-")
sc.pp.calculate_qc_metrics(
    adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True
)
sc.pp.filter_genes(adata, min_cells=3)

adata = adata[adata.obs.pct_counts_mt < 5, :]
adata = adata[adata.obs.n_genes_by_counts < 5000, :]

# Post-QC count matrix (before ref gene intersect / normalize).
# Reloaded later to build CAFs_TotalGene.h5ad; also used by update/precess.py.
os.makedirs(RESULT_DIR, exist_ok=True)
adata.write(f"{RESULT_DIR}/adata_caf.h5ad")
print(f"Wrote post-QC counts: {RESULT_DIR}/adata_caf.h5ad")


# =============================================================================
# 2. Intersect with breast reference; write train / test
# =============================================================================
ref = sc.read_h5ad(REF_H5AD)
over_ = np.intersect1d(ref.var.index, adata.var.index)
ref = ref[:, over_]
adata = adata[:, over_]

sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

# HVGs from concatenated ref+query
test = sc.AnnData.concatenate(ref, adata)
sc.pp.highly_variable_genes(test, n_top_genes=3000)
test = test[:, test.var.highly_variable]

adata_hvg = adata[:, test.var.index]
ref = ref[:, test.var.index]

ref.write("BreastCancer_train.h5ad")
adata_hvg.write("BreastCancer_test.h5ad")


# =============================================================================
# 3.  major annotation (legacy API block)
#    Prefer 4-cellmemory-total.py for the cleaner re-run that writes cls.h5ad
#    (used downstream as cls_total.h5ad).
# =============================================================================
import CellMemory as cellmemory
import numpy as np
import pandas as pd
import scanpy as sc
import os
import sklearn.metrics as metrics
import matplotlib.pyplot as plt

# --- Train (celltype_major) ---
batch_size = 140
Dataset = "data/BreastCancer_train.h5ad"
Label = "celltype_major"
Project = "BreastCancer_major"
cellmemory.train(Dataset, Label=Label, Project=Project, batch_size=batch_size)

# --- Predict (celltype_major) ---
Label = "celltype_major"
Project = "BreastCancer_major"
Model = Project + "_ckpt.pth"
batch_size = 300

test_data = sc.read_h5ad(f"{RESULT_DIR}/BreastCancer_test.h5ad")

os.chdir(RESULT_DIR)
adata_pred, total_attn, AttentionGene = cellmemory.predict(
    test_data,
    Model,
    Project=Project,
    batch_size=batch_size,
    out_Tag=True,
    data_parallel=True,
)

os.chdir(Project)

sc.tl.pca(adata_pred)
sc.pp.neighbors(adata_pred)
sc.tl.umap(adata_pred)
# KNOWN: original used sc.tl.leiden(adata) with undefined `adata` — see KNOWN_ISSUES.md
# sc.tl.leiden(adata)

sc.pl.umap(adata_pred, color="Pred", frameon=False, show=False)
plt.title("Annotation", fontsize=20)
plt.savefig("cls_Pred.pdf", dpi=300, bbox_inches="tight")

sc.pl.umap(adata_pred, color="Pred_prob", frameon=False, show=False)
plt.title("Confidence Score", fontsize=20)
plt.savefig("cls_Pred_prob.pdf", dpi=300, bbox_inches="tight")

adata_pred.obs["sample"] = test_data.obs["sample"]
sc.pl.umap(adata_pred, color="sample", frameon=False, show=False)
plt.title("Sample", fontsize=20)
plt.savefig("cls_sample.pdf", dpi=300, bbox_inches="tight")

adata_pred.write("cls_total.h5ad")

# =============================================================================
# 4. Keep CAF cells only → CAFs_TotalGene.h5ad
# =============================================================================
# Reload post-QC counts from adata_caf.h5ad (written after section 1), then subset
# to CAF barcodes and normalize for total-gene CAF matrix.
# KNOWN: original indexes adata_pred with undefined `pos` after defining pos_caf.
# Logic kept as-is; see KNOWN_ISSUES.md. Intended index is almost certainly pos_caf.
pos_caf = np.where(adata_pred.obs["Pred"] == "CAFs")[0]
adata_pred = adata_pred[pos, :]  # noqa: F821  — original bug; intended: pos_caf
adata = sc.read_h5ad(f"{RESULT_DIR}/adata_caf.h5ad")
adata = adata[pos_caf, :]
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

adata.write("CAFs_TotalGene.h5ad")


# =============================================================================
# 5. EXPLORATORY — supervised CAF-subtype CellMemory (NOT the final paper path)
#    Final CAF subtypes come from unsupervised clustering in 5-caf-cluster.py
#    (Pred_prob > 0.7 + Leiden → m/i/t/d/vCAF).
# =============================================================================
# adata = sc.read_h5ad(
#     f"{RESULT_DIR}/BreastCancer_major/CAFs_TotalGene.h5ad"
# )
# ref = sc.read_h5ad(
#     f"{RESULT_DIR}/BreastCancer_major/CAFs/ref_CAFs_rawCounts.h5ad"
# )
# sc.pp.normalize_total(ref, target_sum=1e4)
# sc.pp.log1p(ref)
# over_ = np.intersect1d(ref.var.index, adata.var.index)
# ref = ref[:, over_]
# adata = adata[:, over_]
# sc.pp.highly_variable_genes(ref, n_top_genes=4000)
# ref = ref[:, ref.var.highly_variable]
# adata = adata[:, ref.var.index]
# ref.write("CAFs_train_hvg4k.h5ad")
# adata.write("CAFs_test_hvg4k.h5ad")
# ... CellMemory train/predict on Label='CAFs_type' ...
