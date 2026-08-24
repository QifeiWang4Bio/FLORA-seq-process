#!/usr/bin/env python3
# =============================================================================
# Module : 5-matrix / Species-mixing QC (human–mouse UMI collision)
# Dataset: FLORA-seq mixR3 (hg38 + mm39 dual-genome Cell Ranger matrix)
# Source : 1_Cell_line/0-bin/6-speices-mix.py
#
# Purpose:
#   Estimate sequencing / multiplet quality from per-cell human vs mouse UMI
#   fractions. Cells with minority-species UMI > 20% are labeled "Mixed".
#
# Expected counts (after min_genes=1000 filter):
#   Human 1140 | Mouse 759 | Mixed 89 (4.48%)
#
# Outputs (under mixR3/outs/):
#   update_fig/Mix_species.png
#   update_fig/Mix_species_.png
#   update_result/species_mixing_statistics.csv
#   update_result/adata_mix.h5ad   # input for 2-cluster.py
# =============================================================================

import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc

WORKDIR = (
    "/path/to/floraseq/cell_line/"
    "mixR3/outs"
)
MIN_GENES = 1000
MIX_FRAC_CUTOFF = 0.20  # minority-species UMI fraction → Mixed


def _qc_by_genome(adata):
    """Split dual-genome matrix and compute per-species QC metrics."""
    human = adata[:, adata.var["genome"] == "hg38"].copy()
    mouse = adata[:, adata.var["genome"] == "mm39"].copy()
    for ad in (human, mouse):
        ad.var["mt"] = ad.var_names.str.startswith("MT-")
        sc.pp.calculate_qc_metrics(
            ad, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True
        )
    return human, mouse


def classify_species(total_h, total_m, mix_cutoff=MIX_FRAC_CUTOFF):
    """
    Assign Human / Mouse / Mixed from per-cell human and mouse UMI totals.

    Rule:
      - majority species by UMI count
      - if minority / (human + mouse) > mix_cutoff → Mixed
    """
    total = total_h + total_m
    # avoid divide-by-zero; empty cells stay Human/Mouse by majority (0 vs 0 → Mouse)
    minority_frac = np.where(
        total_h > total_m,
        np.divide(total_m, total, out=np.zeros(len(total)), where=total > 0),
        np.divide(total_h, total, out=np.zeros(len(total)), where=total > 0),
    )
    labels = np.where(total_h > total_m, "Human", "Mouse").astype(object)
    labels[minority_frac > mix_cutoff] = "Mixed"
    return labels


def print_stats(adata):
    """Print cell counts and per-species UMI / gene summaries."""
    print("\n" + "=" * 60)
    print("SPECIES MIXING STATISTICS")
    print("=" * 60)

    counts = adata.obs["species"].value_counts()
    n = len(adata)
    print(f"\nTotal cells: {n}")
    for sp, c in counts.items():
        print(f"  {sp}: {c} ({100 * c / n:.2f}%)")
    if "Mixed" in counts:
        print(f"  Doublet (Mixed) rate: {100 * counts['Mixed'] / n:.2f}%")

    print("\nPer-species median UMIs / genes:")
    for sp in ["Human", "Mouse", "Mixed"]:
        if sp not in adata.obs["species"].values:
            continue
        mask = adata.obs["species"] == sp
        if sp == "Human":
            umi = adata.obs.loc[mask, "total_counts_human"]
            genes = adata.obs.loc[mask, "n_genes_human"]
        elif sp == "Mouse":
            umi = adata.obs.loc[mask, "total_counts_mouse"]
            genes = adata.obs.loc[mask, "n_genes_mouse"]
        else:
            umi = (
                adata.obs.loc[mask, "total_counts_human"]
                + adata.obs.loc[mask, "total_counts_mouse"]
            )
            genes = (
                adata.obs.loc[mask, "n_genes_human"]
                + adata.obs.loc[mask, "n_genes_mouse"]
            )
        print(
            f"  {sp}: UMI median={umi.median():.0f}, "
            f"genes median={genes.median():.0f}"
        )
    print("=" * 60)


def main():
    os.chdir(WORKDIR)
    os.makedirs("update_fig", exist_ok=True)
    os.makedirs("update_result", exist_ok=True)

    # ----- load dual-genome matrix -----
    data = sc.read_10x_h5("filtered_feature_bc_matrix.h5")
    data.var_names_make_unique()

    # ----- filter cells with too few genes in either genome -----
    data_h, data_m = _qc_by_genome(data)
    sc.pp.filter_cells(data_h, min_genes=MIN_GENES)
    sc.pp.filter_cells(data_m, min_genes=MIN_GENES)
    keep = np.union1d(data_h.obs_names, data_m.obs_names)
    data = data[keep, :].copy()

    # ----- recompute species-wise counts on filtered cells -----
    data_h, data_m = _qc_by_genome(data)
    data.obs["total_counts_human"] = 0.0
    data.obs["total_counts_mouse"] = 0.0
    data.obs["n_genes_human"] = 0.0
    data.obs["n_genes_mouse"] = 0.0
    data.obs.loc[data_h.obs_names, "total_counts_human"] = data_h.obs["total_counts"].values
    data.obs.loc[data_m.obs_names, "total_counts_mouse"] = data_m.obs["total_counts"].values
    data.obs.loc[data_h.obs_names, "n_genes_human"] = data_h.obs["n_genes_by_counts"].values
    data.obs.loc[data_m.obs_names, "n_genes_mouse"] = data_m.obs["n_genes_by_counts"].values

    data.obs["species"] = classify_species(
        data.obs["total_counts_human"].to_numpy(),
        data.obs["total_counts_mouse"].to_numpy(),
    )
    data.obsm["Mix_species"] = np.vstack(
        [data.obs["total_counts_human"], data.obs["total_counts_mouse"]]
    ).T

    # ----- scatter: human UMI vs mouse UMI -----
    data.uns["species_colors"] = ["#1f77b4", "#D3D3D3", "#2ca02c"]
    sc.pl.embedding(data, basis="Mix_species", color="species", show=False)
    plt.savefig("update_fig/Mix_species.png", dpi=300, bbox_inches="tight")
    plt.close()

    # padded axes for publication-style axis limits (dummy cell at corner)
    data_pad = data.concatenate(data[:1].copy())
    data_pad.obs.loc[data_pad.obs_names[-1], "total_counts_human"] = 60000
    data_pad.obs.loc[data_pad.obs_names[-1], "total_counts_mouse"] = 50000
    data_pad.obsm["Mix_species"] = np.vstack(
        [data_pad.obs["total_counts_human"], data_pad.obs["total_counts_mouse"]]
    ).T
    data_pad.uns["species_colors"] = ["#1f77b4", "#D3D3D3", "#2ca02c"]
    sc.pl.embedding(
        data_pad, basis="Mix_species", color="species", show=False, size=60
    )
    plt.savefig("update_fig/Mix_species_.png", dpi=600, bbox_inches="tight")
    plt.close()

    print_stats(data)
    pd.DataFrame(
        {
            "Species": data.obs["species"].value_counts().index,
            "Count": data.obs["species"].value_counts().values,
            "Percentage": 100
            * data.obs["species"].value_counts().values
            / len(data),
        }
    ).to_csv("update_result/species_mixing_statistics.csv", index=False)

    data.write_h5ad("update_result/adata_mix.h5ad")
    print("Wrote update_result/adata_mix.h5ad")


if __name__ == "__main__":
    main()
