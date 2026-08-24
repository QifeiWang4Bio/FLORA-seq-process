#!/usr/bin/env python3
# =============================================================================
# Module : 4-sj / Plot junction-saturation curves (SJ diversity)
# Source : Python block in 1_Cell_line/0-bin/13-transcipt-diverisity.py
#
# Input  : RSeQC *.junctionSaturation_plot.r files under WORK_DIR
# Output : Benchmark_Junctions_{All,Known,Novel}.pdf
#
# Final plot methods (VASA omitted intentionally, matching published figures):
#   10x Genomics - 3', Smart-seq3, FLORA-seq (3'), FLORA-seq (Random)
# =============================================================================

import os
import re

import matplotlib.pyplot as plt
import seaborn as sns

plt.rcParams["pdf.fonttype"] = 42

WORK_DIR = "/path/to/floraseq/benchmark/sj_diversity"

# Color / style: 10x as baseline; Smart-seq3 as gold standard;
# FLORA 3' polyA; FLORA Random highlighted.
FILES = {
    "10x Genomics - 3": {
        "path": "10x_3_human_saturation.junctionSaturation_plot.r",
        "color": "#7f8c8d",
        "marker": "o",
        "style": "--",
    },
    "Smart-seq3": {
        "path": "smartseq3_saturation.junctionSaturation_plot.r",
        "color": "#2ecc71",
        "marker": "^",
        "style": "--",
    },
    "FLORA-seq (3)": {
        "path": "flora_human_3_saturation.junctionSaturation_plot.r",
        "color": "#3498db",
        "marker": "s",
        "style": "-",
    },
    "FLORA-seq (Random)": {
        "path": "flora_human_R_saturation.junctionSaturation_plot.r",
        "color": "#e74c3c",
        "marker": "D",
        "style": "-",
    },
}

PLOTS_CONFIG = [
    ("z", "All Splicing Junctions (Sensitivity)", "Total Junctions", "All"),
    ("y", "Known Splicing Junctions (Annotation)", "Known Junctions", "Known"),
    ("w", "Novel Splicing Junctions (Discovery)", "Novel Junctions", "Novel"),
]


def parse_rseqc_r_file(file_path):
    """Parse RSeQC junctionSaturation_plot.r vectors x/y/z/w."""
    data = {}
    if not os.path.exists(file_path):
        print(f"Warning: File not found - {file_path}")
        return None

    with open(file_path, "r") as f:
        content = f.read()

    patterns = {
        "x": r"x=c\((.*?)\)",  # depth %
        "y": r"y=c\((.*?)\)",  # known
        "z": r"z=c\((.*?)\)",  # all
        "w": r"w=c\((.*?)\)",  # novel
    }
    for key, pattern in patterns.items():
        match = re.search(pattern, content)
        if match:
            data[key] = [int(n) for n in match.group(1).split(",")]
    return data


def plot_all_metrics():
    sns.set_style("ticks")
    sns.set_context("paper", font_scale=1.4)

    for metric_key, title, ylabel, suffix in PLOTS_CONFIG:
        plt.figure(figsize=(7, 6))
        print(f"Plotting {title}...")

        for label, config in FILES.items():
            data = parse_rseqc_r_file(os.path.join(WORK_DIR, config["path"]))
            if not data or metric_key not in data:
                continue

            x = data["x"]
            y = data[metric_key]
            plt.plot(
                x,
                y,
                label=label,
                color=config["color"],
                marker=config["marker"],
                linestyle=config["style"],
                linewidth=2.5,
                markersize=6,
                alpha=0.9,
            )
            print(f"  [{suffix}] {label}: {y[-1]:,}")

        plt.title(title, fontsize=16, fontweight="bold", pad=15)
        plt.xlabel("Resampling Depth (% of Total Reads)", fontsize=14)
        plt.ylabel(f"Number of {ylabel}", fontsize=14)
        plt.grid(True, which="major", linestyle="--", alpha=0.4)
        sns.despine()
        plt.legend(frameon=False, fontsize=12)
        plt.tight_layout()

        save_name = os.path.join(WORK_DIR, f"Benchmark_Junctions_{suffix}.pdf")
        plt.savefig(save_name, dpi=300)
        print(f"Saved: {save_name}\n")
        plt.close()


if __name__ == "__main__":
    os.makedirs(WORK_DIR, exist_ok=True)
    plot_all_metrics()
