#!/usr/bin/env python3
# =============================================================================
# Module : 4-sj / Validate FLORA-seq splice junctions vs long-read SJs
# Dataset: FLORA-seq K562 SJs vs ENCODE K562 long-read SJs
# Source : Python block in 1_Cell_line/0-bin/14-long-read-benchmark.py
#
# Prerequisite: run 3-long-read-benchmark.sh to produce *_sjs.bed files.
#
# Matching rule:
#   A FLORA SJ (support >= MIN_SCORE) is validated if any long-read SJ on the
#   same chrom has intron start/end within TOLERANCE bp (fuzzy match).
#   Long-read reference = union of ENCODE BED files below.
#
# Final figure: Fig_LongRead_Validation_Panel.pdf
# =============================================================================

import os
from collections import defaultdict

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

plt.rcParams["pdf.fonttype"] = 42

# ================= Configuration =================
# Default = 3'+Random (mixR3). Uncomment another FLORA_BED to compare designs.
FLORA_BED = (
    "/path/to/floraseq/cell_line/"
    "mixR3_human/outs/flora_k562_sjs.bed"
)
# FLORA_BED = (
#     "/path/to/floraseq/cell_line/"
#     "mixR_human/outs/flora_k562_sjs.bed"
# )
# FLORA_BED = (
#     "/path/to/floraseq/cell_line/"
#     "mix3_human/outs/flora_k562_sjs.bed"
# )

# Long-read SJ BEDs (union used as reference)
LONG_BED_FILES = [
    "/path/to/floraseq/benchmark/long_read/ENCFF479SQR_sjs.bed",
    "/path/to/floraseq/benchmark/long_read/ENCFF322UJU_sjs.bed",
    "/path/to/floraseq/benchmark/long_read/ENCFF504GVG_sjs.bed",
    "/path/to/floraseq/benchmark/long_read/ENCFF645UVN_sjs.bed",
    "/path/to/floraseq/benchmark/long_read/ENCFF661OEY_sjs.bed",
]
# Optional SGNex-only reference (not used for the final figure):
# LONG_BED_FILES = [
#     "/path/to/floraseq/benchmark/long_read/"
#     "SGNex_K562_PacBio-SMRTcell_replicate7_run1_sjs.bed"
# ]

TOLERANCE = 10  # +/- bp fuzzy match on intron boundaries
MIN_SCORE = 2   # minimum FLORA supporting reads
OUT_DIR = "/path/to/floraseq/benchmark/long_read/fig"
# ================================================


def parse_regtools_bed(bed_path):
    """Parse regtools BED12; return {(chrom, intron_start, intron_end): score}."""
    sjs = defaultdict(int)
    if not os.path.exists(bed_path):
        print(f"Warning: File not found {bed_path}")
        return sjs

    with open(bed_path, "r") as f:
        for line in f:
            if not line.strip():
                continue
            cols = line.strip().split("\t")
            chrom = cols[0]
            start = int(cols[1])
            score = int(cols[4])

            try:
                block_sizes = [int(x) for x in cols[10].strip(",").split(",")]
                block_starts = [int(x) for x in cols[11].strip(",").split(",")]
            except (IndexError, ValueError):
                continue

            intron_start = start + block_sizes[0]
            intron_end = start + block_starts[1]
            sjs[(chrom, intron_start, intron_end)] += score

    return sjs


def bin_support(x):
    if x == 2:
        return "2"
    if x == 3:
        return "3"
    if 4 <= x <= 10:
        return "4-10"
    if x > 10:
        return ">10"
    return "Others"


def main():
    os.makedirs(OUT_DIR, exist_ok=True)
    os.chdir(OUT_DIR)

    print(f"Loading FLORA SJs from {os.path.basename(FLORA_BED)}...")
    flora_sjs = parse_regtools_bed(FLORA_BED)
    flora_sjs = {k: v for k, v in flora_sjs.items() if v >= MIN_SCORE}
    print(f"  - {len(flora_sjs)} high-confidence FLORA SJs (Support >= {MIN_SCORE}).")

    long_index = defaultdict(list)
    total_long_sjs_count = 0

    print("Loading Long-Read reference files...")
    for fpath in LONG_BED_FILES:
        fname = os.path.basename(fpath)
        print(f"  - Processing {fname}...")
        current_sjs = parse_regtools_bed(fpath)
        for chrom, start, end in current_sjs.keys():
            long_index[chrom].append((start, end))
            total_long_sjs_count += 1

    print(f"  - Total Long-Read SJ records loaded: {total_long_sjs_count}")
    for chrom in long_index:
        long_index[chrom] = list(set(long_index[chrom]))

    print("Starting validation (fuzzy matching)...")
    validated_count = 0
    results = []
    total_flora = len(flora_sjs)
    process_counter = 0

    for (chrom, start, end), score in flora_sjs.items():
        process_counter += 1
        if process_counter % 5000 == 0:
            print(f"  - Processed {process_counter}/{total_flora} SJs...")

        is_validated = False
        if chrom in long_index:
            for l_start, l_end in long_index[chrom]:
                if abs(start - l_start) <= TOLERANCE and abs(end - l_end) <= TOLERANCE:
                    is_validated = True
                    break

        if is_validated:
            validated_count += 1
        results.append({"Support": score, "Validated": is_validated})

    val_rate = validated_count / total_flora * 100 if total_flora else 0.0
    print("\n=== Final Result ===")
    print(f"Validated: {validated_count} / {total_flora} ({val_rate:.2f}%)")

    df = pd.DataFrame(results)
    df.to_csv("sj_validation_results_final.csv", index=False)
    print("Results saved to sj_validation_results_final.csv")

    # ---------------- plots ----------------
    print("Generating plots...")
    sns.set_style("whitegrid")
    plt.rcParams.update({"font.size": 12})

    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    ax1 = axes[0]
    n_val = int(df["Validated"].sum())
    n_not = len(df) - n_val
    ax1.pie(
        [n_val, n_not],
        labels=[f"Validated\n({n_val})", f"Not Validated\n({n_not})"],
        autopct="%1.1f%%",
        colors=["#2ecc71", "#e74c3c"],
        startangle=90,
        explode=(0.05, 0),
        textprops={"fontsize": 14},
    )
    ax1.set_title(
        "Global Validation Rate of FLORA-seq SJs\n(vs. ENCODE Long-Reads)",
        fontsize=16,
    )

    ax2 = axes[1]
    df["Support_Bin"] = df["Support"].apply(bin_support)
    bin_order = ["2", "3", "4-10", ">10"]
    summary = df.groupby("Support_Bin")["Validated"].mean().reset_index()
    summary["Validated_Pct"] = summary["Validated"] * 100

    sns.barplot(
        data=summary,
        x="Support_Bin",
        y="Validated_Pct",
        order=bin_order,
        palette="Blues_d",
        ax=ax2,
    )

    summary_ordered = summary.set_index("Support_Bin").reindex(bin_order)
    for i, row in enumerate(summary_ordered.itertuples()):
        if pd.isna(row.Validated_Pct):
            continue
        ax2.text(
            i,
            row.Validated_Pct + 1,
            f"{row.Validated_Pct:.1f}%",
            ha="center",
            color="black",
            fontweight="bold",
        )

    ax2.set_xlabel("Number of FLORA-seq Reads Supporting the SJ", fontsize=14)
    ax2.set_ylabel("Percentage Validated by Long-Reads (%)", fontsize=14)
    ax2.set_title("Higher Support Correlates with Long-Read Validation", fontsize=16)
    ax2.set_ylim(0, 110)

    plt.tight_layout()
    plt.savefig("Fig_LongRead_Validation_Panel.pdf", dpi=300)
    print("Plot saved to Fig_LongRead_Validation_Panel.pdf")


if __name__ == "__main__":
    main()
