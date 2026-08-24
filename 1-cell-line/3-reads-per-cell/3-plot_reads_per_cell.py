#!/usr/bin/env python3
# =============================================================================
# Module : 3-reads-per-cell / Plot saturation curves (final figure)
# Source : plotting block in 1_Cell_line/0-bin/4-reads-per-cell-human-only.sh
#
# Inputs (from 2-reads-per-cell.sh):
#   flora_stats.csv, 10x_3_stats.csv, 10x_5_stats.csv
#   columns: Method, TargetDepth, CellID, UMI_Count, Gene_Count
#
# Outputs:
#   benchmark_genes.pdf  — median genes vs mean reads/cell
#   benchmark_umis.pdf   — median UMIs vs mean reads/cell
# =============================================================================

import os

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

FILES = {
    "FLORA-seq": (
        "/path/to/floraseq/cell_line/"
        "mixR3_human/outs/depth/flora_stats.csv"
    ),
    "10x_3": (
        "/path/to/floraseq/benchmark/10x_3/"
        "293t_mix_1k_v3-1/293t_mix_1k_v3-1_human/outs/depth/10x_3_stats.csv"
    ),
    "10x_5": (
        "/path/to/floraseq/benchmark/10x_5/"
        "10k_hgmm_5pv2_human/outs/depth/10x_5_stats.csv"
    ),
}

OUTDIR = (
    "/path/to/floraseq/cell_line/"
    "mixR3_human/outs/depth"
)


def main():
    dfs = []
    for method, fpath in FILES.items():
        if not os.path.exists(fpath):
            print(f"Missing: {fpath}")
            continue
        dfs.append(pd.read_csv(fpath))

    if not dfs:
        raise SystemExit("No stats CSVs found. Run 2-reads-per-cell.sh first.")

    os.chdir(OUTDIR)
    df = pd.concat(dfs, ignore_index=True)
    df_summary = (
        df.groupby(["Method", "TargetDepth"]).median(numeric_only=True).reset_index()
    )

    sns.set_style("ticks")

    plt.figure(figsize=(6, 5))
    sns.lineplot(
        data=df_summary,
        x="TargetDepth",
        y="Gene_Count",
        hue="Method",
        marker="o",
        linewidth=2,
    )
    plt.xlabel("Mean Reads per Cell")
    plt.ylabel("Median Genes Detected")
    plt.title("Sequencing Saturation (Gene Detection)")
    plt.grid(True, linestyle="--", alpha=0.6)
    plt.legend(frameon=False)
    plt.tight_layout()
    plt.savefig("benchmark_genes.pdf")
    plt.close()

    plt.figure(figsize=(6, 5))
    sns.lineplot(
        data=df_summary,
        x="TargetDepth",
        y="UMI_Count",
        hue="Method",
        marker="o",
        linewidth=2,
    )
    plt.xlabel("Mean Reads per Cell")
    plt.ylabel("Median UMIs Detected")
    plt.title("Sequencing Saturation (Sensitivity)")
    plt.grid(True, linestyle="--", alpha=0.6)
    plt.legend(frameon=False)
    plt.tight_layout()
    plt.savefig("benchmark_umis.pdf")
    plt.close()

    print(f"Wrote {OUTDIR}/benchmark_genes.pdf and benchmark_umis.pdf")


if __name__ == "__main__":
    main()
