#!/usr/bin/env python3
# =============================================================================
# Module : 3-mHypo / count.h5ad + CellMemory train/test HVG matrices
# Source : 3_mHypo/bin/old/ref.py
#
# Prerequisite: 1-process-qc.R → counts.mtx
# FINAL HVG: hvg5k (C7/C25) → 4-anno-coarse; hvg8k (C66) → 9-anno-c66
# Embedded CellMemory train/predict block below is exploratory.
#
# Environment: conda activate sc-py; scanpy
# =============================================================================

import os
import random
import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad

# 2. Query data
os.chdir('/path/to/floraseq/mhypo/seurat_out')
adata = sc.read('counts.mtx')

barcode = pd.read_csv('barcode.csv')
gene = pd.read_csv('features.csv')
adata.obs.index = list(barcode['x'])
adata.var.index = list(gene['x'])
adata.var_names_make_unique()
adata.obs_names_make_unique()

adata.write(filename='count.h5ad')  # Write the count object.
# Filter selected genes.
sc.pp.filter_genes(adata, min_cells=3)

adata.var['mt'] = adata.var_names.str.startswith('MT-')
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

# 1. reference data
# celltype: C25_named(25)   C66_named(66)  
data =sc.read_h5ad('/path/to/external-data/gw/data/mouse_hypothalamic/local.h5ad')
data.var.index = list(data.var['feature_name'])
data.raw = None

data_train = data.copy()
random.seed(123)
random_pos = np.arange(data_train.shape[0])
np.random.shuffle(random_pos)
data_train = data_train[random_pos, :]

over_ = np.intersect1d(data_train.var.index, adata.var.index)
data_train = data_train[:, over_]
adata = adata[:, over_]

sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
adata_ = adata.copy()


# FINAL: hvg5k → C7/C25 (4-anno-coarse)
data_train = data.copy()
random.seed(123)
random_pos = np.arange(data_train.shape[0])
np.random.shuffle(random_pos)
data_train = data_train[random_pos, :]

adata = adata_.copy()
data_train = data_train[:, over_]
adata = adata[:, over_]

sc.pp.highly_variable_genes(data_train, n_top_genes=5000)
data_train = data_train[:, data_train.var.highly_variable]

adata = adata[:, data_train.var.index]

data_train.write('toCellMemory/mHypo_train_hvg5k.h5ad')
adata.write('toCellMemory/mHypo_test_hvg5k.h5ad')

# FINAL: hvg8k → C66 (9-anno-c66)
data_train = data.copy()
random.seed(123)
random_pos = np.arange(data_train.shape[0])
np.random.shuffle(random_pos)
data_train = data_train[random_pos, :]

# adata = adata_.copy()
data_train = data_train[:, over_]
adata = adata[:, over_]

sc.pp.highly_variable_genes(data_train, n_top_genes=8000)
data_train = data_train[:, data_train.var.highly_variable]

adata = adata[:, data_train.var.index]

data_train.write('toCellMemory/mHypo_train_hvg8k.h5ad')
adata.write('toCellMemory/mHypo_test_hvg8k.h5ad')



# EXPLORATORY: CellMemory train/predict sweeps + Pred summary
# (C7/C25 Full_bin training was originally TrainModel.py; C66 final is 9-anno-c66)
import CellMemory as cellmemory
import numpy as np
import pandas as pd
import scanpy as sc
import os
import sklearn.metrics as metrics
import matplotlib.pyplot as plt

# ***************************** #
# >>>          Train        <<< #
# ***************************** #
# hvg5k-C7
batch_size = 60
Dataset = 'data/mHypo_train_hvg5k.h5ad'  # train data path
Label = 'C7_named'   # Supertype  Subclass
Project = 'mHypo_meng_hvg5k_C7'  # Project name
cellmemory.train(Dataset, Label=Label, Project=Project, batch_size=batch_size)

# hvg5k-C25
batch_size = 60
Dataset = 'data/mHypo_train_hvg5k.h5ad'  # train data path
Label = 'C25_named'   # Supertype  Subclass
Project = 'mHypo_meng_hvg5k_C25'  # Project name
cellmemory.train(Dataset, Label=Label, Project=Project, batch_size=batch_size)



# ***************************** #
# >>>        Predict        <<< #
# ***************************** #
import CellMemory as cellmemory
import numpy as np
import pandas as pd
import scanpy as sc
import os
import sklearn.metrics as metrics
import matplotlib.pyplot as plt
import anndata as ad
from scipy.sparse import csr_matrix

group = ['mHypo_meng_hvg5k_C7','mHypo_meng_hvg5k_C25']  # ,'mHypo_meng_hvg5k_C66'
for g in group:
    Project = g # Project name
    Model = Project+'_ckpt.pth'
    batch_size = 150
    # 
    os.chdir('/path/to/floraseq/mhypo/seurat_out')
    test_data = sc.read_h5ad('/path/to/floraseq/mhypo/seurat_out/toCellMemory/mHypo_test_hvg5k.h5ad')
    # 
    adata_pred, total_attn, AttentionGene = \
    cellmemory.predict(test_data, Model, Project=Project, batch_size=batch_size, out_Tag=True, data_parallel=True)
    # 
    os.chdir(Project)
    AttentionGene.to_csv('AttentionGene.csv',index=None)
    adata_attn = ad.AnnData(csr_matrix(np.array(total_attn)))
    adata_attn.obs = test_data.obs
    adata_attn.var = test_data.var
    # adata_attn.write('adata_attn.h5ad')


# Pred
group = ['mHypo_meng_hvg5k_C7','mHypo_meng_hvg5k_C25']

Pred = pd.DataFrame()
Pred.index = list(test_data.obs.index)
for g in group:
    cls = sc.read_h5ad(g+'/cls.h5ad')
    Pred[g] = list(cls.obs['Pred'])

Pred.columns = ['hvg5k_C7','hvg5k_C25']

# Pred prob
group = ['mHypo_meng_hvg5k_C7','mHypo_meng_hvg5k_C25']

Pred_prob = pd.DataFrame()
Pred_prob.index = list(test_data.obs.index)
for g in group:
    cls = sc.read_h5ad(g+'/cls.h5ad')
    Pred_prob[g] = list(cls.obs['Pred_prob'])

Pred_prob.columns = ['hvg5k_C7','hvg5k_C25']

# Overall, the HVG5k setting produced the higher confidence score.
Pred.to_csv('summary_Pred.csv')
Pred_prob.to_csv('summary_Pred_prob.csv')
