#!/usr/bin/env python3
# =============================================================================
# Module : 3-mHypo / Imprinted-gene heatmap at C66 (neurons)
# Source : 3_mHypo/bin/old/update.py
#
# Prerequisite: 10-neuron-c66.py → adata_neuron_norm.h5ad
# Coarse (L1) imprint: 6-imprint-coarse.py
#
# Environment: conda activate sc-py; scanpy
# =============================================================================

import os
import scanpy as sc
import matplotlib.pyplot as plt

plt.rcParams["pdf.fonttype"] = 42
sc.set_figure_params(figsize=(5, 5), frameon=False)

os.chdir(
    "/path/to/floraseq/mhypo/seurat_out/mHypo_meng_hvg8k_filter_C66_bin"
)

adata = sc.read_h5ad("adata_neuron_norm.h5ad")
groupby = "C66" if "C66" in adata.obs.columns else "anno_l3"

marker = [
    "Ube3a", "Mirg", "Asb4", "Th", "Cdkn1c", "Calcr", "Grb10", "Dlk1",
    "Peg3", "Peg10", "Peg12", "Ndn", "Magel2", "Nnat", "Inpp5f",
]
marker = [g for g in marker if g in adata.var_names]

sc.pl.matrixplot(
    adata, marker, groupby,
    colorbar_title="Expression \nscaled by column",
    standard_scale="var", vmin=-1, vmax=1.5, cmap="RdBu_r",
)
plt.savefig("heatmap_hypo_C66_imprint.pdf", dpi=600, bbox_inches="tight")
