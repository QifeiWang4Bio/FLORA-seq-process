#!/usr/bin/env python3
# =============================================================================
# Module : 2-CAF / CellMemory major-type train + predict (cleaner re-run)
# Dataset: Breast cancer FLORA-seq (BreastCancer_train / test from 3-prepare-annotate.py)
# Source : 2_hBreast/bin/old/CellMemory_total.py
#
# Purpose:
#   Re-train / predict celltype_major with CellMemory (bin_scale), write cls.h5ad
#   with Pred / Pred_prob. Downstream scripts refer to this as cls_total.h5ad.
#
# Prerequisite:
#   BreastCancer_train.h5ad
#   BreastCancer_test.h5ad
#
# Section 0 = primary run used for annotation confidence scores
# Section 1 = 5-fold validation variants (exploratory / robustness)
# =============================================================================

import os

import CellMemory as cm
import matplotlib.pyplot as plt
import scanpy as sc

# Paths used in the original re-run (adjust if train/test live elsewhere)
TRAIN_H5AD = "/path/to/data/meng_breast/BreastCancer_train.h5ad"
TEST_H5AD = "/path/to/data/meng_breast/BreastCancer_test.h5ad"

Label = "celltype_major"


# =============================================================================
# 0. Primary run: hvg3k-noFilter-bin-1GPU
#    Output: <Project>/cls.h5ad  (+ Pred / Pred_prob UMAPs)
# =============================================================================
batch_size = 140
Dataset = TRAIN_H5AD
Project = "meng_breast_total_noF_bin"
cm.train(
    Dataset,
    Label=Label,
    Project=Project,
    batch_size=batch_size,
    exp_bin="bin_scale",
)

batch_size = 70
test_data = sc.read_h5ad(TEST_H5AD)
adata_pred = cm.predict(
    test_data,
    Project=Project,
    batch_size=batch_size,
    out_Tag=False,
    exp_bin="bin_scale",
)

sc.pp.pca(adata_pred)
sc.pp.neighbors(adata_pred)
sc.tl.umap(adata_pred)
os.chdir(Project)
sc.set_figure_params(figsize=(5, 5))

sc.pl.umap(
    adata_pred, color="Pred", frameon=False, show=False, use_raw=False, size=10
)
plt.title("CellMemory", fontsize=20)
plt.savefig("cls_Pred.png", dpi=300, bbox_inches="tight")

sc.pl.umap(
    adata_pred,
    color="Pred_prob",
    frameon=False,
    show=False,
    use_raw=False,
    size=10,
)
plt.title("Confidence Score", fontsize=20)
plt.savefig("cls_Pred_prob.png", dpi=300, bbox_inches="tight")

adata_pred.write("cls.h5ad")
# Copy / rename to cls_total.h5ad under BreastCancer_major/ for 5-caf-cluster.py

