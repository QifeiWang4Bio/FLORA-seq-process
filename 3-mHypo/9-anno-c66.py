#!/usr/bin/env python3
# =============================================================================
# Module : 3-mHypo / CellMemory train + predict at C66
# Source : 3_mHypo/bin/old/model_C66.py
# Prerequisite: 2-process-hvg.py → mHypo_train/test_hvg8k.h5ad
# FINAL: mHypo_meng_hvg8k_filter_C66_bin → coembed.h5ad (→ 10-neuron-c66)
# Section 1 (no-bin) is exploratory.
# Environment: CellMemory, scanpy, torch, GPU recommended
# =============================================================================

# CellMemory
import CellMemory as cellmemory
import numpy as np
import pandas as pd
import scanpy as sc
import os
import sklearn.metrics as metrics
import matplotlib.pyplot as plt

def get_attention_list(
                adata_cls, 
                attn_mat, 
                ref_idx2celltype, 
                ref_gene,
                topN=50
):
    # get top attention score gene list
    gene_list = pd.DataFrame()
    for celltype in ref_idx2celltype:
        attn = torch.mean(attn_mat[np.where((adata_cls.obs["Pred"]==celltype))[0],:], dim=0)
        if torch.isnan(attn[0]):
            continue
        attn_sort = torch.sort(attn, descending=True)
        atten_gene = ref_gene[attn_sort[1][0:topN]]
        gene_list[celltype] = list(atten_gene)
    return gene_list

# ---------------------------------------------------------------------------
# EXPLORATORY: filter-noBin (no bins=). Final path is section 2 (filter-Bin).
# ---------------------------------------------------------------------------
# 1) filter-noBin
# ***************************** #
# >>>          Train        <<< #
# ***************************** #
# hvg8k
batch_size = 128
Dataset = '/path/to/data/meng_mHypo/mHypo_train_hvg8k.h5ad'  # train data path
Label = 'C66_named'   # Supertype  Subclass
Project = 'mHypo_meng_hvg8k_filter_C66'  # Project name

cellmemory.train(Dataset, Label=Label, Project=Project, batch_size=batch_size)


# ***************************** #
# >>>        Predict        <<< #
# ***************************** #
Model = Project+'_ckpt.pth'
batch_size = 256
# 
test = sc.read_h5ad('/path/to/data/meng_mHypo/mHypo_test_hvg8k.h5ad')
train = sc.read_h5ad(Dataset)
train.obs['batch_id'] = 'Reference'
test.obs['batch_id'] = 'Query'
test_data = sc.AnnData.concatenate(train, test, index_unique=None)

# 1. test set
adata_pred, total_attn, out_gene_ids =\
cellmemory.predict(test, Model, Project=Project, batch_size=batch_size, out_Tag=True)
# 
import torch
import anndata as ad
from scipy.sparse import csr_matrix
os.chdir(Project)

total_attn = total_attn.numpy()
out_gene_ids = out_gene_ids.numpy()
def make_attn_mat(test_data, total_attn, out_gene_ids):
	attn_mat = np.array(np.full((test_data.shape[0], test_data.shape[1]), 0), dtype=np.float32)
	for i in range(test_data.shape[0]):
		temp = (out_gene_ids[i]-1)
		temp = np.array(temp[~(temp==3000)], dtype=int)
		attn_mat[i, temp] = total_attn[i, 0:len(temp)]
		adata_attn = ad.AnnData(attn_mat)
		adata_attn.obs.index = list(test_data.obs.index)
		adata_attn.var.index = list(test_data.var.index)
		adata_attn.X = csr_matrix(adata_attn.X)
	return adata_attn

adata_attn = make_attn_mat(test_data, total_attn, out_gene_ids)
adata_attn.write('adata_attn.h5ad')

ref_idx2celltype = list(adata_pred.obs['Pred'].value_counts().index)
ref_gene = np.array(pd.read_csv(Project+'/'+Project+"_ref_gene.txt", header=None)[0])
AttentionGene = get_attention_list(adata_pred, total_attn, ref_idx2celltype, ref_gene, topN=100)
AttentionGene.to_csv('AttentionGene.csv',index=None)


# 2. merge set
os.chdir('..')
adata_pred = \
cellmemory.predict(test_data, Model, Project=Project, batch_size=batch_size, out_Tag=False)

sc.set_figure_params(figsize=(5, 5))
sc.pp.pca(adata_pred)
sc.pp.neighbors(adata_pred)
sc.tl.umap(adata_pred)

sc.pl.umap(adata_pred, color='Pred', frameon=False,show=False,use_raw=False, size=1)
plt.title('CellMemory', fontsize=20)
plt.savefig('coembed_pred.png', dpi=300, bbox_inches="tight")

adata_pred.obs['batch_id'] = list(test_data.obs['batch_id'])
sc.pl.umap(adata_pred, color='batch_id', frameon=False,show=False,use_raw=False, size=1)
plt.title('Batch', fontsize=20)
plt.savefig('coembed_batch.png', dpi=300, bbox_inches="tight")

adata_pred.write('coembed.h5ad')



# ---------------------------------------------------------------------------
# FINAL: filter-Bin → Project mHypo_meng_hvg8k_filter_C66_bin / coembed.h5ad
# Consumed by 10-neuron-c66.py
# ---------------------------------------------------------------------------
# 2) filter-Bin
# ***************************** #
# >>>          Train        <<< #
# ***************************** #
# hvg8k
batch_size = 128
Dataset = '/path/to/data/meng_mHypo/mHypo_train_hvg8k.h5ad'  # train data path
Label = 'C66_named'   # Supertype  Subclass
Project = 'mHypo_meng_hvg8k_filter_C66_bin'  # Project name

cellmemory.train(Dataset, Label=Label, Project=Project, batch_size=batch_size, bins='bins')


# ***************************** #
# >>>        Predict        <<< #
# ***************************** #
Model = Project+'_ckpt.pth'
batch_size = 256
# 
test = sc.read_h5ad('/path/to/data/meng_mHypo/mHypo_test_hvg8k.h5ad')
train = sc.read_h5ad(Dataset)
train.obs['batch_id'] = 'Reference'
test.obs['batch_id'] = 'Query'
test_data = sc.AnnData.concatenate(train, test, index_unique=None)

# 1. test set
adata_pred, total_attn, AttentionGene = \
cellmemory.predict(test, Model, Project=Project, batch_size=batch_size, out_Tag=True, bins='bins')
# 
import torch
import anndata as ad
from scipy.sparse import csr_matrix
os.chdir(Project)

total_attn = total_attn.numpy()
out_gene_ids = out_gene_ids.numpy()
def make_attn_mat(test_data, total_attn, out_gene_ids):
	attn_mat = np.array(np.full((test_data.shape[0], test_data.shape[1]), 0), dtype=np.float32)
	for i in range(test_data.shape[0]):
		temp = (out_gene_ids[i]-1)
		temp = np.array(temp[~(temp==3000)], dtype=int)
		attn_mat[i, temp] = total_attn[i, 0:len(temp)]
		adata_attn = ad.AnnData(attn_mat)
		adata_attn.obs.index = list(test_data.obs.index)
		adata_attn.var.index = list(test_data.var.index)
		adata_attn.X = csr_matrix(adata_attn.X)
	return adata_attn

adata_attn = make_attn_mat(test_data, total_attn, out_gene_ids)
adata_attn.write('adata_attn.h5ad')

ref_idx2celltype = list(adata_pred.obs['Pred'].value_counts().index)
ref_gene = np.array(pd.read_csv(Project+'/'+Project+"_ref_gene.txt", header=None)[0])
AttentionGene = get_attention_list(adata_pred, total_attn, ref_idx2celltype, ref_gene, topN=100)
AttentionGene.to_csv('AttentionGene.csv',index=None)

# 2. merge set
os.chdir('..')
adata_pred = \
cellmemory.predict(test_data, Model, Project=Project, batch_size=batch_size, out_Tag=False, bins='bins')

sc.set_figure_params(figsize=(5, 5))
sc.pp.pca(adata_pred)
sc.pp.neighbors(adata_pred)
sc.tl.umap(adata_pred)

sc.pl.umap(adata_pred, color='Pred', frameon=False,show=False,use_raw=False, size=1)
plt.title('CellMemory', fontsize=20)
plt.savefig('coembed_pred.png', dpi=300, bbox_inches="tight")

sc.pl.umap(adata_pred, color='batch_id', frameon=False,show=False,use_raw=False, size=1)
plt.title('Batch', fontsize=20)
plt.savefig('coembed_batch.png', dpi=300, bbox_inches="tight")

adata_pred.write('coembed.h5ad')











