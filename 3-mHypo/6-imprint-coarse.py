#!/usr/bin/env python3
# =============================================================================
# Module : 3-mHypo / Coarse marker UMAPs + imprinted genes (L1)
# Source : marker_plot.py; actions from 7-fig7.py §1 / §3.1
# Prerequisite: 4-anno-coarse.py → cls_run.h5ad / exp_run.h5ad
# C66 imprint: 11-imprint-c66.py
# Environment: conda activate sc-py; scanpy
# =============================================================================

import os
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt

try:
    from cellmemory.plot import createCorlormap
except ImportError:
    createCorlormap = None

sc.set_figure_params(figsize=(5, 5), frameon=False)
plt.rcParams["pdf.fonttype"] = 42

C25_DIR = "/path/to/floraseq/mhypo/seurat_out/mHypo_meng_hvg5k_C25"
C7_DIR = "/path/to/floraseq/mhypo/seurat_out/mHypo_meng_hvg5k_C7"
REF_H5AD = "/path/to/external-data/gw/data/mouse_hypothalamic/local.h5ad"
IMPRINT_GENE = os.path.join(C7_DIR, "Imprinted_gene.txt")
IMPRINT_FROM = os.path.join(C7_DIR, "Imprinted_gene_from.txt")


def _cmap(name):
    return createCorlormap(name) if createCorlormap is not None else ("viridis" if name == "exp" else "magma")


def _rotate_xticks(fig):
    plt.sca(plt.get(fig, "axes")["mainplot_ax"])
    plt.xticks(
        rotation=45,
        rotation_mode="anchor",
        verticalalignment="center",
        horizontalalignment="right",
    )


# --- coarse marker UMAPs ---
os.chdir(C7_DIR)
os.makedirs("marker", exist_ok=True)

cls = sc.read_h5ad(os.path.join(C25_DIR, "cls_run.h5ad"))
adata = sc.read_h5ad(os.path.join(C25_DIR, "exp_run.h5ad"))  # already norm + coarse
adata.obsm["X_umap"] = cls.obsm["X_umap"]

gs = [
    "Snap25", "Syt1", "Slc17a6", "Gad2", "Pdgfra", "Cspg4", "C1ql1", "Vcan",
    "Mbp", "Plp1", "Mobp", "Ptgds", "Aqp4", "Agt", "Gja1", "Atp1b2", "Ccdc153",
    "Rax", "Cdhr4", "Scn7a", "Slco1a4", "Vtn", "Slc6a13", "Igfbp7", "Luc7l3",
]
for gene in gs:
    if gene not in adata.var_names:
        continue
    sc.pl.umap(
        adata, color=gene, frameon=False, show=False, size=7,
        palette="tab10", cmap=_cmap("exp"),
    )
    plt.title(gene, fontsize=20)
    plt.savefig(f"marker/cls_C7_{gene}.png", dpi=300, bbox_inches="tight")


# --- imprint (L1) ---
os.makedirs(os.path.join(C7_DIR, "imprinted"), exist_ok=True)
os.chdir(os.path.join(C7_DIR, "imprinted"))

# curated panel (Meng list)
im_heatmap = [
    "Axl", "Nlrp2",  # As
    "Pon3", "Ccdc40", "Nnat",  # Ep
    "Th", "Htr2a", "Kcnk9", "Dact2", "Zdbf2", "Snrpn", "Dlk1", "Peg10",  # Neu
    "Ascl2", "Mkrn3",  # OPC
    "Gatm", "Galnt6",  # Oli
    "Dcn",  # Vas
]
im_heatmap = [g for g in im_heatmap if g in adata.var_names]

sc.pl.matrixplot(
    adata, im_heatmap, "coarse",
    colorbar_title="Expression \nscaled by column",
    standard_scale="var", vmin=-1, vmax=1.5, cmap="RdBu_r",
)
plt.savefig("heatmap_bulk_imprint.pdf", dpi=600, bbox_inches="tight")

im_umap = ["Axl", "Galnt6", "Kcnk9", "Peg10"]
for gene in im_umap:
    if gene not in adata.var_names:
        continue
    sc.pl.umap(adata, color=gene, frameon=False, show=False, size=10, cmap=_cmap("attn"))
    plt.title(gene, fontsize=20)
    plt.savefig(f"cls_imprint_select_{gene}.png", dpi=300, bbox_inches="tight")

# paternal / maternal from gene lists
gs_imprint = pd.read_csv(IMPRINT_GENE, header=None)
gs_imprint_from = pd.read_csv(IMPRINT_FROM, header=None)
gs_imprint.columns = ["imprint"]
gs_imprint["imprint_from"] = gs_imprint_from[0]

gs_imprint_paternal = np.intersect1d(
    list(gs_imprint[gs_imprint.imprint_from == "Paternal"].imprint),
    list(adata.var.index),
)
gs_imprint_maternal = np.intersect1d(
    list(gs_imprint[gs_imprint.imprint_from == "Maternal"].imprint),
    list(adata.var.index),
)

sc.pl.matrixplot(
    adata, gs_imprint_paternal, "coarse",
    colorbar_title="Expression \nscaled by column",
    standard_scale="var", vmin=-1, vmax=1.5, cmap="RdBu_r",
)
plt.savefig("heatmap_bulk_imprint_paternal.pdf", dpi=600, bbox_inches="tight")

sc.pl.matrixplot(
    adata, gs_imprint_maternal, "coarse",
    colorbar_title="Expression \nscaled by column",
    standard_scale="var", vmin=-1, vmax=1.5, cmap="RdBu_r",
)
plt.savefig("heatmap_bulk_imprint_maternal.pdf", dpi=600, bbox_inches="tight")



# --- select imprint panel ---
im_select = [
    "Axl", "Nlrp2", "Cd81",  # As
    "Pon3", "Ccdc40", "Nnat",  # Ep
    "Th", "Htr2a", "Kcnk9", "Mirg", "Zdbf2", "Snrpn", "Dlk1", "Peg10",  # Neu
    "Ampd3", "Slc22a3", "Smoc1",  # OPC
    "Gatm", "Galnt6", "Gab1",  # Oli
    "Dcn",  # Vas
]
im_select = [g for g in im_select if g in adata.var_names]

sc.pl.matrixplot(
    adata, im_select, "coarse",
    colorbar_title="Expression \nscaled by column",
    standard_scale="var", vmin=-1, vmax=1.5, cmap="RdBu_r",
)
plt.savefig("heatmap_bulk_imprint_select.pdf", dpi=600, bbox_inches="tight")
