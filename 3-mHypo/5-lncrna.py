#!/usr/bin/env python3
# =============================================================================
# Module : 3-mHypo / lncRNA heatmap at coarse level
# Source : 3_mHypo/bin/old/lncRNA.py
# Prerequisite: 4-anno-coarse.py → exp_run.h5ad
# Environment: conda activate sc-py; scanpy
# =============================================================================

import os
import random
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt

# coarse-level
lnc = pd.read_csv('/path/to/reference/ensembl/mouse/lncRNA_gene.txt',header=None)
os.chdir('/path/to/floraseq/mhypo/seurat_out/mHypo_meng_hvg5k_C25')
adata = sc.read_h5ad('exp_run.h5ad')
os.makedirs('lncRNA', exist_ok=True)
os.chdir('lncRNA')

lnc_gene = np.intersect1d(list(adata.var.index), lnc[0].astype(str))
adata_lnc = adata[:, lnc_gene]
adata_lnc.uns['log1p']['base'] = None
sc.pp.highly_variable_genes(adata_lnc, n_top_genes=3000)
adata_lnc_hvg = adata_lnc[:, adata_lnc.var.highly_variable]

sc.tl.rank_genes_groups(adata_lnc_hvg, 'coarse', method='wilcoxon', tie_correct=True, n_genes=20)  # , rankby_abs=True

gene = []
gene_pd = pd.DataFrame()
gene_dict = {}
for i in ['Astrocytes', 'Ependymal-like', 'Neurons', 'OPC', 'Oligodendrocytes', 'Vascular']:
    temp = list(adata_lnc_hvg.uns['rank_genes_groups']['names'][i][0:12])
    gene.extend(temp)
    gene_pd[i] = temp
    gene_dict[i] = temp

# Retain at most 500 cells per cell type for the single-cell heatmap.
group = list(adata_lnc.obs.coarse.value_counts().index)
pos = []
for g in group:
    temp = np.where(adata_lnc.obs['coarse']==g)[0]
    if len(temp) > 500:
        temp = list(np.random.choice(temp, 500))
    pos.extend(temp)

adata_scale_filter = adata_lnc[pos, :]

# single cell
sc.pl.heatmap(adata_scale_filter, gene, groupby='coarse', cmap='viridis', dendrogram=False,standard_scale='var',swap_axes=True,\
    vmin=0, vmax=1)
plt.savefig('heatmap_sc_coarse_lnRNA_test.pdf', dpi=300, bbox_inches="tight")

gene_pd.to_csv('lncRNA_coarse_top12.csv',index=None)
