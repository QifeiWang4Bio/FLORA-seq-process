#!/usr/bin/env python3
# =============================================================================
# Module : 5-matrix / Extract HEK293T matrices across platforms
# Source : 1_Cell_line/0-bin/9-process-matrix.py
#
# Purpose:
#   After FLORA-seq cell-type annotation (2-cluster.py), identify HEK293T
#   barcodes for each benchmarking platform and subset human-mapped count
#   matrices. Outputs feed a fair reads-per-cell / depth benchmark.
#
# Platforms:
#   1) FLORA-seq  — HEK293T from adata_celltype.h5ad
#   2) 10x-3'     — Leiden annotation on dual-genome matrix
#   3) 10x-5'     — Leiden annotation on dual-genome matrix
#   4) VASA-seq   — HEK293T plate counts
#   5) SMART-seq3 — zUMIs expression matrix
#
# Outputs (under .../benchmark/depth_293t/):
#   floraseq_293t.h5ad, 10x_3_293t.h5ad, 10x_5_293t.h5ad,
#   vasaseq_293t.h5ad, smartseq3_293t.h5ad
#
# Also writes hek293t_barcode.txt next to each dual-genome Cell Ranger run
# (used later by subset-bam in 3-reads-per-cell/).
# =============================================================================

import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp

HEK293T_MARKERS = [
    "hg38_VIM",
    "hg38_ANXA1",
    "hg38_ANXA2",
    "hg38_LGALS3",
    "hg38_ITGA6",
    "hg38_AMOT",
]
OUTDIR = "/path/to/floraseq/benchmark/depth_293t"


def print_qc(adata, name):
    adata = adata.copy()
    adata.var["mt"] = adata.var_names.str.startswith("MT-")
    sc.pp.calculate_qc_metrics(
        adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True
    )
    print(
        f"{name}: cells={adata.n_obs}, "
        f"median genes={np.median(adata.obs['n_genes_by_counts']):.1f}, "
        f"median UMIs={np.median(adata.obs['total_counts']):.1f}"
    )
    return adata


def cluster_dual_genome(workdir, hek_leiden, nih_leiden, fig_prefix="fig"):
    """
    Load dual-genome 10x matrix, cluster, annotate HEK293T vs NIH3T3,
    write hek293t_barcode.txt, return barcode Series.
    """
    os.chdir(workdir)
    os.makedirs(fig_prefix, exist_ok=True)

    data = sc.read_10x_h5("filtered_feature_bc_matrix.h5")
    data.var_names_make_unique()
    data.var["mt"] = data.var_names.str.startswith("MT-")
    sc.pp.calculate_qc_metrics(
        data, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True
    )
    print(
        f"[{workdir}] all cells={data.n_obs}, "
        f"median genes={np.median(data.obs['n_genes_by_counts']):.1f}, "
        f"median UMIs={np.median(data.obs['total_counts']):.1f}"
    )

    sc.pp.normalize_total(data, target_sum=1e6)
    sc.pp.log1p(data)
    sc.pp.highly_variable_genes(data, n_top_genes=2000)
    data_hvg = data[:, data.var.highly_variable].copy()
    sc.pp.pca(data_hvg)
    sc.pp.neighbors(data_hvg)
    sc.tl.umap(data_hvg)
    sc.tl.leiden(data_hvg)

    sc.pl.umap(data_hvg, color="leiden", frameon=False, show=False, size=30)
    plt.title("Leiden Clustering", fontsize=14)
    plt.savefig(f"{fig_prefix}/umap_leiden.png", dpi=300, bbox_inches="tight")
    plt.close()

    data.obsm["X_umap"] = data_hvg.obsm["X_umap"]
    sc.pl.umap(data, color=HEK293T_MARKERS, frameon=False, show=False, size=62, ncols=5)
    plt.savefig(f"{fig_prefix}/umap_markers.png", dpi=300, bbox_inches="tight")
    plt.close()

    celltype = np.array(["Unknown"] * data_hvg.n_obs, dtype=object)
    celltype[data_hvg.obs["leiden"].isin(hek_leiden)] = "HEK293T"
    celltype[data_hvg.obs["leiden"].isin(nih_leiden)] = "NIH3T3"
    data_hvg.obs["celltype"] = celltype

    hek = data_hvg[data_hvg.obs["celltype"] == "HEK293T"]
    pd.DataFrame(hek.obs_names).to_csv("hek293t_barcode.txt", index=False, header=False)
    print(f"  Wrote hek293t_barcode.txt ({hek.n_obs} cells)")
    return pd.read_csv("hek293t_barcode.txt", header=None)[0]


def subset_human_matrix(human_outdir, barcodes):
    """Subset a human-only Cell Ranger matrix by barcode whitelist."""
    os.chdir(human_outdir)
    data = sc.read_10x_h5("filtered_feature_bc_matrix.h5")
    data.var_names_make_unique()
    keep = np.intersect1d(barcodes.astype(str), data.obs_names.astype(str))
    return data[keep, :].copy()


# =============================================================================
# 1. FLORA-seq — HEK293T from annotated mixR3 object
# =============================================================================
flora_mix = "/path/to/floraseq/cell_line/mixR3/outs"
os.chdir(flora_mix)

adata_ct = sc.read_h5ad("update_result/adata_celltype.h5ad")
hek = adata_ct[adata_ct.obs["celltype"] == "HEK293T"].copy()
# keep -1 suffix for Cell Ranger matrix matching
pd.DataFrame(hek.obs_names).to_csv("hek293t_barcode.txt", index=False, header=False)
print(f"FLORA-seq HEK293T barcodes: {hek.n_obs}")

bar_flora = pd.read_csv("hek293t_barcode.txt", header=None)[0]
data_1 = subset_human_matrix(
    "/path/to/floraseq/cell_line/mixR3_human/outs",
    bar_flora,
)


# =============================================================================
# 2. 10x-3' — cluster dual-genome, then subset human matrix
# =============================================================================
bar_10x3 = cluster_dual_genome(
    workdir=(
        "/path/to/floraseq/benchmark/10x_3/"
        "293t_mix_1k_v3-1/293t_mix_1k_v3-1_run/outs"
    ),
    hek_leiden=["0", "1", "5"],
    nih_leiden=["2", "3", "4", "6", "7"],
)
data_2 = subset_human_matrix(
    (
        "/path/to/floraseq/benchmark/10x_3/"
        "293t_mix_1k_v3-1/293t_mix_1k_v3-1_human/outs"
    ),
    bar_10x3,
)


# =============================================================================
# 3. 10x-5' — cluster dual-genome, then subset human matrix
# =============================================================================
bar_10x5 = cluster_dual_genome(
    workdir=(
        "/path/to/floraseq/benchmark/10x_5/"
        "10k_hgmm_5pv2_run/outs"
    ),
    hek_leiden=["0", "1", "6", "8", "9", "11"],
    nih_leiden=["2", "3", "4", "5", "7", "10", "12"],
)
data_3 = subset_human_matrix(
    "/path/to/floraseq/benchmark/10x_5/10k_hgmm_5pv2_human/outs",
    bar_10x5,
)


# =============================================================================
# 4. VASA-seq — HEK293T plate transcript counts
# =============================================================================
os.chdir(
    "/path/to/floraseq/benchmark/vasaseq/"
    "vasaseq_output_HEK293T/06_counts/"
)
adata = sc.read_csv(
    "SRR14783059_S1_L001_uniaggGenes_total.TranscriptCounts.tsv",
    delimiter="\t",
    first_column_names=True,
).T

idx = adata.var.index.to_series()
adata.var["gene_id"] = idx.str.extract(r"^([^_]+)", expand=False)
adata.var["gene_name"] = idx.str.extract(r"^[^_]+_([^_]+)", expand=False)

adata = adata[:, adata.var["gene_name"].notna()].copy()
gene_names = adata.var["gene_name"].values
unique_genes, inverse_idx = np.unique(gene_names, return_inverse=True)
X = adata.X if sp.issparse(adata.X) else sp.csr_matrix(adata.X)
agg_mat = sp.csr_matrix(
    (np.ones(len(inverse_idx)), (np.arange(len(inverse_idx)), inverse_idx)),
    shape=(len(inverse_idx), len(unique_genes)),
)
data_4 = sc.AnnData(
    X=X @ agg_mat,
    obs=adata.obs.copy(),
    var=pd.DataFrame(index=unique_genes),
)
data_4.var["gene_name"] = data_4.var.index
sc.pp.calculate_qc_metrics(data_4, inplace=True)
print(
    f"VASA-seq: cells={data_4.n_obs}, "
    f"median genes={np.median(data_4.obs['n_genes_by_counts']):.1f}, "
    f"median UMIs={np.median(data_4.obs['total_counts']):.1f}"
)


# =============================================================================
# 5. SMART-seq3 — zUMIs expression export
# =============================================================================
# Prerequisite: export_for_python.R under smartseq3 zUMIs output
data_dir = (
    "/path/to/floraseq/benchmark/smartseq3/"
    "zumi_fwdprimer/zUMIs_output/expression/for_python"
)
adata = sc.read_mtx(os.path.join(data_dir, "matrix.mtx")).T
genes = pd.read_csv(
    os.path.join(data_dir, "genes.tsv"),
    sep="\t",
    header=None,
    names=["gene_id", "gene_name"],
)
barcodes = pd.read_csv(
    os.path.join(data_dir, "barcodes.tsv"),
    sep="\t",
    header=None,
    names=["barcode"],
)
adata.var_names = genes["gene_id"].values
adata.var["gene_name"] = genes["gene_name"].values
adata.obs_names = barcodes["barcode"].values
sc.pp.calculate_qc_metrics(adata, inplace=True)
print(
    f"SMART-seq3: cells={adata.n_obs}, "
    f"median genes={adata.obs['n_genes_by_counts'].median():.0f}, "
    f"mean UMIs={adata.obs['total_counts'].mean():.0f}"
)
adata.write(os.path.join(data_dir, "umi_matrix.h5ad"))
data_5 = adata.copy()


# =============================================================================
# 6. QC summary + write HEK293T matrices
# =============================================================================
os.makedirs(OUTDIR, exist_ok=True)
data_1 = print_qc(data_1, "FLORA-seq HEK293T")
data_2 = print_qc(data_2, "10x-3 HEK293T")
data_3 = print_qc(data_3, "10x-5 HEK293T")
data_4 = print_qc(data_4, "VASA-seq HEK293T")
data_5 = print_qc(data_5, "SMART-seq3 HEK293T")

os.chdir(OUTDIR)
data_1.write("floraseq_293t.h5ad")
data_2.write("10x_3_293t.h5ad")
data_3.write("10x_5_293t.h5ad")
data_4.write("vasaseq_293t.h5ad")
data_5.write("smartseq3_293t.h5ad")
print(f"Wrote HEK293T matrices to {OUTDIR}")
