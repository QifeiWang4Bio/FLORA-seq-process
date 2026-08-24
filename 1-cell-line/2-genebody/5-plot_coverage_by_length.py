#!/usr/bin/env python3
# =============================================================================

import matplotlib.pyplot as plt
import pandas as pd

cover_1 = pd.read_csv(
    "/path/to/floraseq/benchmark/coverage/"
    "RseQC_floraseq_3k_long.geneBodyCoverage.txt",
    sep="\t",
    index_col=0,
)
cover_2 = pd.read_csv(
    "/path/to/floraseq/benchmark/coverage/"
    "RseQC_floraseq_3k_short.geneBodyCoverage.txt",
    sep="\t",
    index_col=0,
)

cover_1.index = ["FLORA-seq-3-long", "FLORA-seq-R-long"]
cover_2.index = ["FLORA-seq-3-short", "FLORA-seq-R-short"]
cover_merge = pd.concat([cover_1, cover_2])
cover_merge.to_csv(
    "/path/to/floraseq/benchmark/coverage/"
    "genebody_coverage_floraseq_3k_len.csv",
    sep="\t",
)

cover_prop = cover_merge.div(cover_merge.sum(axis=1), axis=0)
plt.rcParams["pdf.fonttype"] = 42

plt.figure(figsize=(10, 6))
for tech in cover_prop.index:
    # highlight FLORA library designs; length strata keep higher alpha
    alpha_val = 1.0 if tech in ["FLORA-seq-3", "FLORA-seq-R"] else 0.8
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
    "genebody_coverage_floraseq_3k_len.pdf",
    dpi=300,
    bbox_inches="tight",
)
plt.close()
