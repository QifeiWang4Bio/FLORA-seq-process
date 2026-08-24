#!/usr/bin/env python3
# ============================================================================

import matplotlib.pyplot as plt
import pandas as pd

cover_merge = pd.read_csv(
    "/path/to/floraseq/benchmark/coverage/"
    "RseQC_5tools.geneBodyCoverage.txt",
    sep="\t",
    index_col=0,
)

# Must match BAM order in cover_5tools.txt
cover_merge.index = [
    "FLORA-seq-3",
    "FLORA-seq-R",
    "VASA-seq",
    "10x-3",
    "10x-5",
    "Smartseq3.5prime_UMI",
    "Smartseq3.Internal",
]

# Final benchmark set
cover_merge = cover_merge.drop(index=["10x-3", "10x-5", "Smartseq3.5prime_UMI"])

cover_prop = cover_merge.div(cover_merge.sum(axis=1), axis=0)
plt.rcParams["pdf.fonttype"] = 42

plt.figure(figsize=(10, 6))
for tech in cover_prop.index:
    alpha_val = 1.0 if tech in ["FLORA-seq-3", "FLORA-seq-R"] else 0.2
    plt.plot(
        cover_prop.columns,
        cover_prop.loc[tech],
        label=tech,
        linewidth=2,
        alpha=alpha_val,
    )

plt.xlabel("Gene body percentile (%)", fontsize=12)
plt.ylabel("Coverage proportion", fontsize=12)
plt.title("Gene Body Coverage Distribution", fontsize=14)
plt.legend(bbox_to_anchor=(1.05, 1), loc="upper left")
plt.xticks(
    [1, 11, 21, 31, 41, 51, 61, 71, 81, 91, 100],
    [0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100],
)
plt.tight_layout()
plt.savefig(
    "/path/to/floraseq/benchmark/coverage/"
    "genebody_coverage_plot_filter.pdf",
    dpi=300,
    bbox_inches="tight",
)
plt.close()
