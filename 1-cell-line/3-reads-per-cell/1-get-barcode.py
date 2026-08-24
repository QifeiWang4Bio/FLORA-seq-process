#!/usr/bin/env python3
# =============================================================================
# Module : 3-reads-per-cell / Export HEK293T barcodes for subset-bam
# Dataset: FLORA-seq, 10x-3', 10x-5' (final reads-per-cell figure)
#
# Prerequisite:
#   ../5-matrix/2-cluster.py          → adata_celltype.h5ad (FLORA)
#   ../5-matrix/3-process-matrix.py   → hek293t_barcode.txt (10x platforms)
# =============================================================================

import os

import pandas as pd
import scanpy as sc

# ---------------------------------------------------------------------------
# 1. FLORA-seq — HEK293T from annotated mixR3 object
# ---------------------------------------------------------------------------
flora_dir = (
    "/path/to/floraseq/cell_line/mixR3/outs"
)
os.chdir(flora_dir)
adata = sc.read_h5ad("update_result/adata_celltype.h5ad")
hek = adata[adata.obs["celltype"] == "HEK293T"].copy()
pd.DataFrame(hek.obs_names).to_csv("hek293t_barcode.txt", index=False, header=False)
print(f"FLORA-seq HEK293T barcodes: {hek.n_obs} → {flora_dir}/hek293t_barcode.txt")

# ---------------------------------------------------------------------------
# 2. 10x-3' — barcodes from 5-matrix/3-process-matrix.py (Leiden annotation)
# ---------------------------------------------------------------------------
bc_10x3 = (
    "/path/to/floraseq/benchmark/10x_3/"
    "293t_mix_1k_v3-1/293t_mix_1k_v3-1_run/outs/hek293t_barcode.txt"
)
n_10x3 = sum(1 for _ in open(bc_10x3))
print(f"10x-3 HEK293T barcodes: {n_10x3} → {bc_10x3}")

# ---------------------------------------------------------------------------
# 3. 10x-5' — barcodes from 5-matrix/3-process-matrix.py (Leiden annotation)
# ---------------------------------------------------------------------------
bc_10x5 = (
    "/path/to/floraseq/benchmark/10x_5/"
    "10k_hgmm_5pv2_run/outs/hek293t_barcode.txt"
)
n_10x5 = sum(1 for _ in open(bc_10x5))
print(f"10x-5 HEK293T barcodes: {n_10x5} → {bc_10x5}")
