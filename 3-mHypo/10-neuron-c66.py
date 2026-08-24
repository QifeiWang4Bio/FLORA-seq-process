#!/usr/bin/env python3
# =============================================================================
# Module : 3-mHypo / C66 neuron UMAP, markers, GABA/Glu neuropeptide plots
# Source : 3_mHypo/bin/old/run_C66.py
#
# Prerequisite: 9-anno-c66.py → coembed.h5ad
# Outputs: figures; adata_neuron_norm.h5ad (→ 11-imprint / 12-interaction-prep)
#
# Environment: conda activate sc-py; scanpy
# =============================================================================

import os
import random
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt

try:
    from cellmemory.plot import createCorlormap
except ImportError:
    createCorlormap = lambda *_a, **_k: "viridis"

# /path/to/floraseq/mhypo/seurat_out/mHypo_meng_hvg8k_filter_C66_bin
C66_DIR = "/path/to/floraseq/mhypo/seurat_out/mHypo_meng_hvg8k_filter_C66_bin"
COUNT = "/path/to/floraseq/mhypo/seurat_out/count.h5ad"

DROP = {
    "C66-53: Etnppl.Astrocytes", "C66-51: Tanycytes",
    "C66-55: Ermn.Oligodendrocytes", "C66-57: OPC", "C66-52: Ependymal",
    "C66-64: Endothelial", "C66-56: Gpr17.Oligodendrocytes",
    "C66-65: Dcn.Fibroblasts", "C66-63: Mural", "C66-54: Lgals3.Astrocytes",
    "C66-58: Tmem119.Immune", "C66-62: ParsTuber", "UnAssigned",
    "C66-59: F13a1.Immune",
}

os.chdir(C66_DIR)
sc.set_figure_params(figsize=(5, 5))
plt.rcParams["pdf.fonttype"] = 42

adata_pred = sc.read_h5ad("coembed.h5ad")
random.seed(123)
adata_pred = adata_pred[np.random.permutation(adata_pred.shape[0]), :]

sc.pl.umap(adata_pred, color="Pred", frameon=False, show=False, use_raw=False, size=1)
plt.title("CellMemory", fontsize=20)
plt.savefig("coembed_pred.png", dpi=600, bbox_inches="tight")

adata_pred.uns["batch_id_colors"] = ["#ff7f0e", "#D3D3D3"]
sc.pl.umap(adata_pred, color="batch_id", frameon=False, show=False, use_raw=False, size=1.5)
plt.title("Batch", fontsize=20)
plt.savefig("coembed_batch.png", dpi=600, bbox_inches="tight")

cls = adata_pred[adata_pred.obs.batch_id.isin(["Query"]), :].copy()
sc.pl.umap(cls, color="Pred", frameon=False, show=False, use_raw=False, size=4)
plt.title("CellMemory", fontsize=20)
plt.savefig("cls_pred.png", dpi=600, bbox_inches="tight")

sc.pl.umap(
    cls, color="Pred", frameon=False, show=False, use_raw=False, size=4,
    legend_loc="on data", legend_fontsize=2,
)
plt.title("CellMemory", fontsize=20)
plt.savefig("cls_pred_legend.png", dpi=600, bbox_inches="tight")

# normalized expression on Query cells
adata = sc.read_h5ad(COUNT)
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
adata = adata[cls.obs.index, :].copy()
adata.obsm["X_umap"] = cls.obsm["X_umap"]

cls_n = cls[~cls.obs["Pred"].isin(DROP), :].copy()
adata_n = adata[cls_n.obs.index, :].copy()
adata_n.obs["C66"] = cls_n.obs["Pred"]
adata_n.obs["anno_l3"] = cls_n.obs["Pred"]
if "Pred_colors" in cls_n.uns:
    adata_n.uns["C66_colors"] = cls_n.uns["Pred_colors"]
adata_n.obsm["X_umap"] = cls_n.obsm["X_umap"]
adata_n.obs["Pred"] = adata_n.obs["C66"]

# for 11-imprint-c66.py / 12-interaction-prep.py
adata_n.write("adata_neuron_norm.h5ad")

# marker genes parsed from C66 names (skip GLU-/Mixed-/non-neuron tokens)
gs = []
for i in list(cls_n.obs.Pred.cat.categories):
    temp = i.split(" ")[1].split(".")[0]
    if temp.startswith(("GLU", "Mixed", "Ependymal", "Tanycytes", "OPC", "ParsTuber", "Mural", "Endothelial")):
        continue
    gs.append(temp)

sc.pl.matrixplot(
    adata_n, gs, "C66", colorbar_title="Expression \nscaled by column",
    standard_scale="var", vmin=-1, vmax=1.5, cmap="RdBu_r",
)
plt.savefig("heatmap_hypo_C66.png", dpi=600, bbox_inches="tight")

os.makedirs("umap_marker", exist_ok=True)
celltype_ = [
    "C66-7: Foxb1.GLU-2", "C66-25: Foxb1.GLU-8", "C66-6: Samd3.GLU-2",
    "C66-24: Pmch.GLU-7", "C66-39: Lhx6.GABA-1", "C66-35: Lef1.GABA-1",
    "C66-40: Hdc.GABA-1",
]
skip = {"C66-1: GLU-1", "C66-23: Mixed.GLU-6", "C66-28: Mixed.GABA-2", "C66-16: Mixed.GLU-4"}
for celltype in celltype_:
    if celltype in skip:
        continue
    gs_gene = celltype.split(" ")[1].split(".")[0]
    adata_n.obs["C66"] = adata_n.obs["Pred"].astype(object)
    adata_n.obs.loc[adata_n.obs["C66"] != celltype, "C66"] = np.nan
    sc.pl.umap(adata_n, color="C66", frameon=False, show=False, use_raw=False, size=8)
    plt.title(celltype, fontsize=20)
    plt.savefig("umap_marker/cls_" + celltype.replace(": ", "_") + ".png", dpi=600, bbox_inches="tight")
    sc.pl.umap(
        adata_n, color=gs_gene, frameon=False, show=False, use_raw=False, size=8,
        cmap=createCorlormap("exp"),
    )
    plt.title(gs_gene, fontsize=20)
    plt.savefig("umap_marker/cls_exp_" + gs_gene + ".png", dpi=600, bbox_inches="tight")

# restore C66 labels after highlight loop
adata_n.obs["C66"] = adata_n.obs["Pred"]

adata_n_f = adata_n[
    adata_n.obs.C66.isin([
        "C66-17: Lpar1.GLU-4", "C66-18: Rfx4.GLU-4", "C66-19: Pomc.GLU-5",
        "C66-20: Tac2.GLU-5", "C66-21: Bace2.GLU-5", "C66-24: Pmch.GLU-7",
        "C66-44: Nkx2-4.GABA-3", "C66-45: Ghrh.GABA-3", "C66-46: Agrp.GABA-4",
        "C66-47: Sst.GABA-4", "C66-49: Satb2.GABA-6",
    ]),
    :,
].copy()

gs_f = []
for i in list(adata_n_f.obs.C66.cat.categories):
    temp = i.split(" ")[1].split(".")[0]
    if temp.startswith(("GLU", "Mixed", "Ependymal", "Tanycytes", "OPC", "ParsTuber", "Mural", "Endothelial")):
        continue
    gs_f.append(temp)

sc.pl.matrixplot(
    adata_n_f, gs_f, "C66", colorbar_title="Expression \nscaled by column",
    standard_scale="var", vmin=-1, vmax=1.5, cmap="RdBu_r",
)
plt.savefig("heatmap_hypo_C66_filter.png", dpi=600, bbox_inches="tight")

# neuropeptide panel (N.txt next to coembed)
gs = pd.read_csv(os.path.join(C66_DIR, "N.txt"), sep="\t", header=None)
gs_2 = list(gs[1][0:21])

temp = []
for i in range(adata_n.shape[0]):
    lab = adata_n.obs.C66.iloc[i]
    if lab == "C66-1: GLU-1":
        temp.append("GLU")
    else:
        temp.append(str(lab).split(".")[1].split("-")[0])
adata_n.obs["Neu"] = temp
adata_n_gaba = adata_n[adata_n.obs.Neu.isin(["GABA"]), :]

fig1 = sc.pl.matrixplot(
    adata_n_gaba, gs_2, "C66", colorbar_title="Expression \nscaled by column",
    standard_scale="var", vmin=-1, vmax=1.5, cmap="RdBu_r", return_fig=True,
)
plt.sca(plt.get(fig1, "axes")["mainplot_ax"])
plt.xticks(rotation=45, rotation_mode="anchor", verticalalignment="center", horizontalalignment="right")
plt.savefig("heatmap_hypo_C66_GABA_N.pdf", dpi=600, bbox_inches="tight")

sc.pl.umap(cls_n, color="Pred", frameon=False, show=False, use_raw=False, size=4)
plt.title("CellMemory (C66)", fontsize=20)
plt.savefig("cls_pred_C66_n.png", dpi=600, bbox_inches="tight")

# optional SJ meta at C66 (main SJ track uses C25 meta from 4-anno-coarse)
bar, batch, merge, celltype = [], [], [], []
for i in range(cls_n.obs.shape[0]):
    bar.append(cls_n.obs.index[i].split("-")[0])
    batch.append(cls_n.obs.index[i].split("_")[1])
    merge.append("group" + batch[i] + "_3r_" + bar[i])
    celltype.append(cls_n.obs.Pred.iloc[i].split(": ")[1])
out = pd.DataFrame({"cell_id": merge, "batch": batch, "barcode": bar, "celltype": celltype})
out.index = out["cell_id"]
out.to_csv("meta.csv", sep="\t")
