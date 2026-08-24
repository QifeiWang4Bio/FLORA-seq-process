#!/usr/bin/env python3
# =============================================================================
# Module : 3-mHypo / Prep C66 neurons for NeuronChat (h5ad → MTX)
# Source : 3_mHypo/bin/6-fig6-interaction.R (§1)
#
# Prerequisite: 9-anno-c66.py → coembed.h5ad; preferably 10-neuron-c66.py
#               → adata_neuron_norm.h5ad (else rebuilt here)
# Next: Rscript 13-interaction-c66.R
#
# Environment: conda activate sc-py; scanpy; scipy
# =============================================================================

import os
import pandas as pd
import scanpy as sc
import scipy.io as sio
import scipy.sparse as sp

C66_DIR = "/path/to/floraseq/mhypo/seurat_out/mHypo_meng_hvg8k_filter_C66_bin"
COUNT = "/path/to/floraseq/mhypo/seurat_out/count.h5ad"
OUT = os.path.join(C66_DIR, "interaction")

DROP = {
    "C66-53: Etnppl.Astrocytes",
    "C66-51: Tanycytes",
    "C66-55: Ermn.Oligodendrocytes",
    "C66-57: OPC",
    "C66-52: Ependymal",
    "C66-64: Endothelial",
    "C66-56: Gpr17.Oligodendrocytes",
    "C66-65: Dcn.Fibroblasts",
    "C66-63: Mural",
    "C66-54: Lgals3.Astrocytes",
    "C66-58: Tmem119.Immune",
    "C66-62: ParsTuber",
    "UnAssigned",
    "C66-59: F13a1.Immune",
}

os.makedirs(OUT, exist_ok=True)
os.chdir(C66_DIR)

neu_path = os.path.join(C66_DIR, "adata_neuron_norm.h5ad")
if os.path.exists(neu_path):
    adata_neu = sc.read_h5ad(neu_path)
else:
    adata_pred = sc.read_h5ad("coembed.h5ad")
    cls = (
        adata_pred[adata_pred.obs.batch_id.isin(["Query"]), :].copy()
        if "batch_id" in adata_pred.obs
        else adata_pred.copy()
    )
    cls_neu = cls[~cls.obs["Pred"].isin(DROP), :].copy()
    adata = sc.read_h5ad(COUNT)
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    adata_neu = adata[cls_neu.obs.index, :].copy()
    adata_neu.obs["C66"] = list(cls_neu.obs["Pred"])
    adata_neu.obs["anno_l3"] = list(cls_neu.obs["Pred"])
    if "X_umap" in cls_neu.obsm:
        adata_neu.obsm["X_umap"] = cls_neu.obsm["X_umap"]
    adata_neu.write(neu_path)

if "anno_l3" not in adata_neu.obs.columns:
    if "C66" in adata_neu.obs.columns:
        adata_neu.obs["anno_l3"] = adata_neu.obs["C66"]
    elif "Pred" in adata_neu.obs.columns:
        adata_neu.obs["anno_l3"] = adata_neu.obs["Pred"]
    else:
        raise KeyError("Need anno_l3 / C66 / Pred in obs")

expr_matrix = adata_neu.X.T
if not sp.issparse(expr_matrix):
    expr_matrix = sp.csr_matrix(expr_matrix)
sio.mmwrite(os.path.join(OUT, "matrix.mtx"), expr_matrix)
pd.DataFrame(adata_neu.var_names).to_csv(
    os.path.join(OUT, "genes.tsv"), sep="\t", header=False, index=False
)
pd.DataFrame(adata_neu.obs_names).to_csv(
    os.path.join(OUT, "barcodes.tsv"), sep="\t", header=False, index=False
)
adata_neu.obs.to_csv(os.path.join(OUT, "metadata.csv"))
print(f"Exported {adata_neu.n_obs} cells → {OUT}/")
print("Next: Rscript 13-interaction-c66.R")
