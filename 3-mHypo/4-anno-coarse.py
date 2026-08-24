#!/usr/bin/env python3
# =============================================================================
# Module : 3-mHypo / Join C7+C25 → coarse; extract neurons; write meta
# Source : 3_mHypo/bin/old/cellmemory.py
#
# Prerequisite: CellMemory cls.h5ad at C7 and C25 (see 2-process-hvg.py
#               exploratory train, or original TrainModel Full_bin).
# Not a C25→C7 label collapse — hybrid `coarse` from both Pred layers.
#
# Outputs (mHypo_meng_hvg5k_C25/):
#   cls_run.h5ad, exp_run.h5ad  → 5-lncrna / 6-imprint-coarse
#   meta.csv (+ neuron/meta.csv) → 8-sj-marvel.R
#
# Environment: conda activate sc-py; scanpy
# =============================================================================

import os
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
from pandas.api.types import CategoricalDtype

try:
    from cellmemory.plot import createCorlormap
except ImportError:
    createCorlormap = lambda *_a, **_k: "viridis"

plt.rcParams["pdf.fonttype"] = 42
sc.set_figure_params(figsize=(5, 5))


os.chdir('/path/to/floraseq/mhypo/seurat_out/mHypo_meng_hvg5k_C25')
cls = sc.read_h5ad('cls.h5ad')
test_data = sc.read_h5ad('/path/to/floraseq/mhypo/seurat_out/toCellMemory/mHypo_test_hvg5k.h5ad')

sc.pl.umap(cls, color='Pred', frameon=False, show=False, size=6)
plt.title('CellMemory', fontsize=20)
plt.savefig('cls_Pred_C25.png', dpi=300, bbox_inches="tight")

sc.pl.umap(cls, color='Pred', frameon=False, show=False, legend_loc='on data', size=6, legend_fontsize=5)
plt.title('CellMemory', fontsize=20)
plt.savefig('cls_Pred_C25_legend.png', dpi=300, bbox_inches="tight")

sc.pl.umap(cls, color='Pred_prob', frameon=False, show=False, size=6)
plt.title('Confidence Score', fontsize=20)
plt.savefig('cls_Pred_C25_prob.png', dpi=300, bbox_inches="tight")

# C7 major cell classes.
cls_c7 = sc.read_h5ad('/path/to/floraseq/mhypo/seurat_out/mHypo_meng_hvg5k_C7/cls.h5ad')
cls.obs['C7'] = cls_c7.obs['Pred']
cls.obs['C7_prob'] = cls_c7.obs['Pred_prob']

sc.pl.umap(cls, color='C7', frameon=False, show=False, size=6)
plt.title('CellMemory (C7)', fontsize=20)
plt.savefig('cls_Pred_C7_.png', dpi=600, bbox_inches="tight")

sc.pl.umap(cls, color='C7_prob', frameon=False, show=False, size=6)
plt.title('CellMemory (C7-Confidence Score)', fontsize=20)
plt.savefig('cls_Pred_C7_prob.png', dpi=300, bbox_inches="tight")

# ***** coarse-level
cls.obs['coarse'] = list(cls.obs['C7'])
cls = cls[np.where((cls.obs['C7']!='C7-5: Immune')&(cls.obs['C7']!='C7-6: ParsTuber'))[0], :]
pos_n = np.where((cls.obs['C7']=='C7-1: GLU')|(cls.obs['C7']=='C7-2: GABA'))[0]
pos_a = np.where((cls.obs['Pred']=='C25-18: Astrocytes'))[0]
pos_ep = np.where((cls.obs['Pred']=='C25-17: Ependymal-like'))[0]
pos_ol = np.where((cls.obs['Pred']=='C25-19: Oligodendrocytes'))[0]
pos_opc = np.where((cls.obs['Pred']=='C25-20: OPC'))[0]
pos_v = np.where((cls.obs['C7']=='C7-7: Vascular'))[0]

cls.obs['coarse'].iloc[pos_n] = 'Neurons'
cls.obs['coarse'].iloc[pos_a] = 'Astrocytes'
cls.obs['coarse'].iloc[pos_ep] = 'Ependymal-like'
cls.obs['coarse'].iloc[pos_ol] = 'Oligodendrocytes'
cls.obs['coarse'].iloc[pos_opc] = 'OPC'
cls.obs['coarse'].iloc[pos_v] = 'Vascular'

pos_1 = np.where((cls.obs['coarse']=='C7-4: Oligo+Precursor'))[0]
pos_2 = np.where((cls.obs['coarse']=='C7-3: Astro-Ependymal'))[0]
cls.obs['coarse'].iloc[pos_1] = 'Neurons'
cls.obs['coarse'].iloc[pos_2] = 'Neurons'

sc.pl.umap(cls, color='coarse', frameon=False, show=False, size=7, palette='tab10')
plt.title('Coarse', fontsize=20)
plt.savefig('cls_coarse_.png', dpi=600, bbox_inches="tight")

# total counts
adata = sc.read_h5ad('/path/to/floraseq/mhypo/seurat_out/count.h5ad')
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
adata = adata[cls.obs.index, :]
adata.obsm['X_umap'] = cls.obsm['X_umap']
adata.obs['coarse'] = cls.obs['coarse']

# Downstream 5-lncrna / 6-imprint-coarse expect these on disk
cls.write('cls_run.h5ad')
adata.write('exp_run.h5ad')

sc.tl.rank_genes_groups(adata, 'coarse', method='wilcoxon')
sc.pl.rank_genes_groups(adata, n_genes=10, sharey=False)

# subsample ≤500 cells / type for heatmap
group = list(adata.obs.coarse.value_counts().index)
pos = []
for g in group:
    temp = np.where(adata.obs['coarse'] == g)[0]
    if len(temp) > 500:
        temp = list(np.random.choice(temp, 500))
    pos.extend(temp)

adata_scale_filter = adata[pos, :].copy()
sc.pp.scale(adata_scale_filter)

gene = []
for g in ['Neurons', 'OPC', 'Oligodendrocytes', 'Astrocytes', 'Ependymal-like', 'Vascular']:
    if g in adata.uns['rank_genes_groups']['names'].dtype.names:
        gene.extend(list(adata.uns['rank_genes_groups']['names'][g][:5]))

adata_scale_filter.obs['coarse'] = adata_scale_filter.obs['coarse'].astype(
    CategoricalDtype(
        categories=['Neurons', 'OPC', 'Oligodendrocytes', 'Astrocytes', 'Ependymal-like', 'Vascular'],
        ordered=True,
    )
)
sc.pl.heatmap(
    adata_scale_filter, gene, groupby='coarse', cmap='viridis',
    dendrogram=False, standard_scale='var', swap_axes=True, vmin=0, vmax=1,
)
plt.savefig('heatmap_sc_coarse.pdf', dpi=300, bbox_inches="tight")


# fine-level (C25-neuron)
pos = np.where((cls.obs['Pred']!='C25-18: Astrocytes')&(cls.obs['Pred']!='C25-17: Ependymal-like')&(cls.obs['Pred']!='C25-19: Oligodendrocytes') \
    &(cls.obs['Pred']!='C25-20: OPC')&(cls.obs['Pred']!='C25-24: Mural+Endothelial')&(cls.obs['Pred']!='C25-23: ParsTuber')\
    &(cls.obs['Pred']!='C25-25: Fibroblasts')&(cls.obs['Pred']!='C25-21: Immune'))[0]

cls_n = cls[pos, :]
adata_n = adata[cls_n.obs.index, :]
cls_n.obsm['coembed_cls'] = cls_n.obsm['X_umap']

sc.tl.pca(cls_n)
sc.pp.neighbors(cls_n)
sc.tl.umap(cls_n)

os.makedirs(
    '/path/to/floraseq/mhypo/seurat_out/mHypo_meng_hvg5k_C25/neuron',
    exist_ok=True,
)
os.chdir('/path/to/floraseq/mhypo/seurat_out/mHypo_meng_hvg5k_C25/neuron')

# Short labels for plots (e.g. "C25-11: GABA-1" → "GABA-1"); keep full for meta
cls_n.obs['Pred'] = list(cls_n.obs['Pred'])
pred_full = list(cls_n.obs['Pred'])
group = list(cls_n.obs.Pred.value_counts().index)
for g in group:
    pos = np.where(cls_n.obs['Pred'] == g)[0]
    cls_n.obs['Pred'].iloc[pos] = g.split(' ')[1]


cls_n.uns['Pred_colors'] = ['#1f77b4', '#aec7e8', '#ffbb78', '#2ca02c', '#d62728', '#ff9896', \
'#8c564b','#f7b6d2','#c5b0d5',  '#e377c2', '#c7c7c7', '#bcbd22', '#17becf', '#9edae5']
sc.pl.umap(cls_n, color='Pred', frameon=False, show=False, size=8)
plt.title('CellMemory (fine-level)', fontsize=20)
plt.savefig('cls_neuron_Pred.png', dpi=600, bbox_inches="tight")

# heatmap
adata_n_scale = adata_n.copy()
sc.pp.scale(adata_n_scale)
adata_n_scale.obs['Pred'] = cls_n.obs['Pred']
gene = ['Foxb1','B230323A14Rik','Pitx2','Cpne9','Grm2','Epha8','Fign','Postn','Nts']

fig = sc.pl.stacked_violin(adata, 'Luc7l3', groupby='Pred',cmap='Reds',swap_axes=False, use_raw=False,standard_scale='var',return_fig=True)
plt.sca(plt.get(fig,'axes')['mainplot_ax'])
plt.xticks(rotation=45,rotation_mode='anchor',verticalalignment='center',horizontalalignment='right')
plt.savefig('vln_Luc7l3.pdf', dpi=300, bbox_inches="tight")

sc.pl.matrixplot(adata_n_scale, gene, groupby='Pred', cmap='Blues',standard_scale='var',swap_axes=False,\
 colorbar_title='column scaled\nExpression')  #    dendrogram =True
plt.savefig('heatmap_average_fine.png', dpi=600, bbox_inches="tight")

gene = ['Lhx6','Crh','Lhx8',
'Vipr2','Vip', 'Pmfbp1',
'Ghrh','Gsx1','Slc18a3',
'Agrp','Sst','Npy',
'Bcl11b','Crym','Ankrd63',
'Th','Slc6a3','Slc18a2',
'Ntng1','Shox2','Tcf7l2',
'Trh','Irx5','Samd3',
'Gpr149','Fezf1','Nr5a1',
'Rfx4','Hcrt','Lhx9',
'Pomc','Tac2','Rxfp1',
'Avp','Oxt','Pou3f2',
'Pmch', 'Otx1', 'Mup6', 
'Foxb1', 'Grm2', 'Pitx2',
]
sc.pl.dotplot(adata_n_scale, gene, groupby='Pred', use_raw=False,standard_scale='var',cmap='Blues',)
plt.savefig('dotplot_average_fine.pdf', dpi=300, bbox_inches="tight")

sc.pl.matrixplot(adata_n_scale, gene, groupby='Pred', cmap='Blues',standard_scale='var',swap_axes=False,\
 colorbar_title='column scaled\nExpression')  #    dendrogram =True
plt.savefig('heatmap_average_fine.pdf', dpi=300, bbox_inches="tight")

gene_ = ['Lhx6','Vipr2','Ghrh','Agrp','Bcl11b','Th','Ntng1','Trh','Gpr149','Rfx4','Pomc','Avp','Pmch','Foxb1']
adata_n.obsm['X_umap'] = cls_n.obsm['X_umap']

for g in gene_:    
    sc.pl.umap(adata_n, color=g, frameon=False, show=False, size=7,cmap='viridis_r')
    plt.title(g, fontsize=20)
    plt.savefig('cls_marker_'+g+'.png', dpi=300, bbox_inches="tight")

g='Gad1'
sc.pl.umap(adata_n, color=g, frameon=False, show=False, size=9, cmap=createCorlormap('exp'))
plt.title(g, fontsize=20)
plt.savefig('cls_marker_'+g+'.pdf', dpi=300, bbox_inches="tight")



# MARVEL meta (C25 neuron labels) — 8-sj-marvel.R reads C25/meta.csv
bar, batch, merge, celltype = [], [], [], []
for i in range(cls_n.obs.shape[0]):
    bar.append(cls_n.obs.index[i].split('-')[0])
    batch.append(cls_n.obs.index[i].split('_')[1])
    merge.append('group' + batch[i] + '_3r_' + bar[i])
    lab = pred_full[i]
    celltype.append(lab.split(': ')[1] if ': ' in lab else lab)

out = pd.DataFrame({'cell_id': merge, 'batch': batch, 'barcode': bar, 'celltype': celltype})
out.index = out['cell_id']
out.to_csv('meta.csv', sep='\t')
out.to_csv('../meta.csv', sep='\t')  # path used by 8-sj-marvel.R

