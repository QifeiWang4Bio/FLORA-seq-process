#!/usr/bin/env python3
# =============================================================================
# Module : 5-matrix / Cluster and annotate FLORA-seq cell lines
# Dataset: FLORA-seq mixR3 (from 1-species-mix.py → adata_mix.h5ad)
# Source : 1_Cell_line/0-bin/7-cluster.py
#
# Workflow:
#   1) Normalize / HVG / PCA / neighbors / UMAP / Leiden
#   2) Marker inspection (cell-line gene panels)
#   3) Assign cell types from Leiden clusters
#   4) Save adata_celltype.h5ad (used to export HEK293T barcodes)
# =============================================================================

import os

import matplotlib.pyplot as plt
import numpy as np
import scanpy as sc

WORKDIR = (
    "/path/to/floraseq/cell_line/"
    "mixR3/outs"
)

# Selected markers for UMAP inspection / matrixplot

CELLTYPE_MARKERS = {
    "HEK293T": ["hg38_ANXA2", "hg38_ANXA1", "hg38_LGALS3"],
    "HeLa": ["hg38_TSIX", "hg38_IRS4", "hg38_DACH1"],
    "K562": ["hg38_GATA1", "hg38_GYPA", "hg38_HBA1"],
    "NIH3T3": ["mm39_Col1a1", "mm39_Col1a2", "mm39_Pdgfra"],
}

# Leiden cluster IDs → cell type (from marker-guided annotation)
LEIDEN_TO_CELLTYPE = {
    "NIH3T3": ["0", "4", "5", "6"],
    "K562": ["1", "7"],
    "HeLa": ["3"],
    "HEK293T": ["2"],
}


def assign_celltype(leiden):
    labels = np.array(["Unknown"] * len(leiden), dtype=object)
    for ctype, clusters in LEIDEN_TO_CELLTYPE.items():
        labels[leiden.isin(clusters)] = ctype
    return labels


def main():
    os.chdir(WORKDIR)
    os.makedirs("update_fig", exist_ok=True)
    os.makedirs("update_result", exist_ok=True)
    sc.set_figure_params(figsize=(5, 5), frameon=False)

    # ----- cluster -----
    adata = sc.read_h5ad("update_result/adata_mix.h5ad")
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, n_top_genes=2000)

    adata_hvg = adata[:, adata.var.highly_variable].copy()
    sc.tl.pca(adata_hvg)
    sc.pp.neighbors(adata_hvg)
    sc.tl.umap(adata_hvg)
    sc.tl.leiden(adata_hvg, resolution=0.5)

    sc.pl.umap(adata_hvg, color="leiden", frameon=False, show=False, size=62)
    plt.title("Leiden Clustering", fontsize=14)
    plt.savefig("update_fig/umap_leiden.png", dpi=300, bbox_inches="tight")
    plt.close()

    # DEG table (optional marker discovery)
    sc.tl.rank_genes_groups(adata_hvg, "leiden", method="wilcoxon")
    deg_df = sc.get.rank_genes_groups_df(adata_hvg, group=None)
    mean_exp = adata_hvg.to_df().groupby(adata_hvg.obs["leiden"]).mean().T
    mean_exp.columns = [f"mean_{c}" for c in mean_exp.columns]
    deg_df = deg_df.merge(mean_exp, left_on="names", right_index=True, how="left")
    deg_sig = deg_df[
        (deg_df["pvals_adj"] < 0.05)
        & (deg_df["logfoldchanges"] > 0.25)
        & (deg_df.filter(like="mean_").max(axis=1) > 1)
    ].copy()
    top = (
        deg_sig.sort_values(["group", "logfoldchanges"], ascending=[True, False])
        .groupby("group")
        .head(30)
        .copy()
    )
    top["rank"] = top.groupby("group").cumcount() + 1
    top.pivot(index="rank", columns="group", values="names").to_csv(
        "update_result/marker_hvg2k_leiden_wilcoxon.csv"
    )

    # ----- assign cell types -----
    adata_hvg.obs["celltype"] = assign_celltype(adata_hvg.obs["leiden"])
    adata_hvg.uns["celltype_colors"] = ["#f3ae61", "#b5aad5", "#d882af", "#91ccae"]
    sc.pl.umap(adata_hvg, color="celltype", frameon=False, show=False, size=62)
    plt.title("Cell Type", fontsize=14)
    plt.savefig("update_fig/umap_celltype.png", dpi=300, bbox_inches="tight")
    plt.close()

    adata.obs["celltype"] = adata_hvg.obs["celltype"].copy()
    sc.pl.matrixplot(
        adata,
        CELLTYPE_MARKERS,
        "celltype",
        colorbar_title="Expression \nscaled by column",
        standard_scale="var",
        vmin=-1,
        vmax=1.5,
        cmap="RdBu_r",
    )
    plt.savefig("update_fig/Heatmap_celltype.png", dpi=300, bbox_inches="tight")
    plt.close()

    # save annotation onto raw counts object
    adata_raw = sc.read_h5ad("update_result/adata_mix.h5ad")
    adata_raw.obs["celltype"] = adata.obs["celltype"].copy()
    adata_raw.write("update_result/adata_celltype.h5ad")

    # ----- 4. QC by cell type (raw counts) -----
    adata_raw.var["mt"] = adata_raw.var_names.str.startswith("MT-")
    sc.pp.calculate_qc_metrics(
        adata_raw, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True
    )
    adata_raw.uns["celltype_colors"] = ["#f3ae61", "#b5aad5", "#d882af", "#91ccae"]
    plt.rcParams["pdf.fonttype"] = 42

    sc.pl.violin(
        adata_raw, "total_counts", groupby="celltype", stripplot=False, inner="box"
    )
    plt.gcf().set_size_inches(h=5, w=5)
    plt.grid(False)
    plt.savefig("update_fig/qc_UMICounts_celltype.png", dpi=300, bbox_inches="tight")
    plt.close()

    sc.pl.violin(
        adata_raw,
        "n_genes_by_counts",
        groupby="celltype",
        stripplot=False,
        inner="box",
    )
    plt.gcf().set_size_inches(h=5, w=5)
    plt.grid(False)
    plt.savefig("update_fig/qc_GeneCounts_celltype.png", dpi=300, bbox_inches="tight")
    plt.close()

    print("Median genes / UMIs by cell type:")
    for ctype in adata_raw.obs["celltype"].value_counts().index:
        sub = adata_raw.obs.loc[adata_raw.obs["celltype"] == ctype]
        print(
            f"  {ctype}: Genes={sub['n_genes_by_counts'].median():.1f}, "
            f"UMIs={sub['total_counts'].median():.1f}"
        )


if __name__ == "__main__":
    main()
