#!/usr/bin/env python3
# =============================================================================
# Module : 6-biotype / Per-cell biotype fraction histograms
# Source : embedded Python block in 1_Cell_line/0-bin/8-biotype-per-cell.sh
#
# Final figure used in the paper:
#   biotype_percell_hist_median.pdf
#   (also writes biotype_percell_hist_mean.pdf)
#
# Layout: 2 rows (protein_coding, lncRNA) × 2 cols (FLORA-seq, 10x-3)
# x = per-cell fraction of reads of that biotype; y = number of cells
# =============================================================================

import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np
import pandas as pd

plt.rcParams["pdf.fonttype"] = 42

DATASETS = {
    "FLORA-seq": (
        "/path/to/floraseq/cell_line/"
        "mixR3_human/outs/biotype_percell/flora_read_biotype.csv"
    ),
    "10x-3": (
        "/path/to/floraseq/benchmark/10x_3/"
        "293t_mix_1k_v3-1/293t_mix_1k_v3-1_human/outs/biotype_percell/"
        "10x3_read_biotype.csv"
    ),
}
BIOTYPES = ["protein_coding", "lncRNA"]
COLORS = {"FLORA-seq": "#E87C47", "10x-3": "#4A90D9"}
BINS = np.arange(0, 1.001, 0.001)
OUT_DIR = (
    "/path/to/floraseq/cell_line/"
    "mixR3_human/outs/biotype_percell"
)


def cell_fractions(path: str, biotype: str) -> np.ndarray:
    """Per-cell fraction of reads assigned to `biotype`."""
    df = pd.read_csv(path)
    total = df.groupby("CB").size()
    hits = (
        df[df["biotype"] == biotype]
        .groupby("CB")
        .size()
        .reindex(total.index, fill_value=0)
    )
    return (hits / total).values


def draw_hist(stat: str) -> None:
    """2×2 histogram grid; annotate each panel with median or mean."""
    stat_fn = np.median if stat == "median" else np.mean
    fig, axes = plt.subplots(
        len(BIOTYPES),
        len(DATASETS),
        figsize=(8 * len(DATASETS), 2 * len(BIOTYPES)),
        sharey="row",
        sharex=True,
    )
    for col, (name, path) in enumerate(DATASETS.items()):
        for row, bt in enumerate(BIOTYPES):
            ax = axes[row, col]
            frac = cell_fractions(path, bt)
            ax.hist(frac, bins=BINS, color=COLORS[name], edgecolor="none")
            val = stat_fn(frac)
            ax.axvline(
                val,
                color="black",
                lw=1.2,
                ls="--",
                label=f"{stat} = {val:.3f}",
            )
            ax.set_title(f"{name}  —  {bt}", fontsize=12)
            ax.set_xlabel("Fraction of reads", fontsize=11)
            if col == 0:
                ax.set_ylabel("Number of cells", fontsize=11)
            ax.set_xlim(0, 1)
            ax.xaxis.set_major_locator(mticker.MultipleLocator(0.2))
            ax.xaxis.set_minor_locator(mticker.MultipleLocator(0.05))
            ax.spines[["top", "right"]].set_visible(False)
            ax.legend(fontsize=9, frameon=False)

    plt.tight_layout()
    out = f"{OUT_DIR}/biotype_percell_hist_{stat}.pdf"
    plt.savefig(out, dpi=300, bbox_inches="tight")
    print(f"Saved: {out}")
    plt.close()


if __name__ == "__main__":
    draw_hist("median")  # final figure
    draw_hist("mean")
