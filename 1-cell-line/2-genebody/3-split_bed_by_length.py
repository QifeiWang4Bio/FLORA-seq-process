#!/usr/bin/env python3
# =============================================================================
# Module : 2-genebody / Split RSeQC gene BED by transcript length
# Source : 1_Cell_line/0-bin/10-split_bed_by_length.py
#
# Used by 4-coverage_by_length.sh to generate short/long (or short/medium/long)
# reference BEDs for FLORA-seq genebody coverage stratified by transcript length.
#
# Examples:
#   THRESHOLDS = [3000]        -> lt3000 / ge3000
#   THRESHOLDS = [2000]        -> lt2000 / ge2000
#   THRESHOLDS = [1000, 5000]  -> short / medium / long (with BIN_LABELS)
# =============================================================================

import os

# ================= Configuration =================
INPUT_BED = "/path/to/reference/RSeQC_ref/hg38_GENCODE.v38.bed"
OUT_PREFIX = "hg38_GENCODE.v38"

# Flexible threshold configuration (unit: bp)
# Examples:
#   [3000]            -> 2 bins: <3kb, >=3kb
#   [1000, 5000]      -> 3 bins: <1kb, 1-5kb, >5kb
#   [500, 2000, 5000] -> 4 bins: <500, 500-2k, 2-5k, >5k
THRESHOLDS = [3000]

# Custom bin labels (optional, leave empty for auto-generated labels)
# If provided, length must be len(THRESHOLDS) + 1
BIN_LABELS = []  # e.g., ["short", "medium", "long"] for [1000, 5000]
# =================================================


def generate_bin_info(thresholds):
    """Generate bin labels and descriptions based on thresholds."""
    if not thresholds:
        raise ValueError("THRESHOLDS cannot be empty")

    thresholds = sorted(thresholds)
    bins = []

    # First bin: < threshold[0]
    bins.append({
        "label": f"lt{thresholds[0]}",
        "desc": f"<{thresholds[0]}bp",
        "range": (0, thresholds[0]),
    })

    # Middle bins: threshold[i] to threshold[i+1]
    for i in range(len(thresholds) - 1):
        bins.append({
            "label": f"{thresholds[i]}to{thresholds[i + 1]}",
            "desc": f"{thresholds[i]}-{thresholds[i + 1]}bp",
            "range": (thresholds[i], thresholds[i + 1]),
        })

    # Last bin: >= threshold[-1]
    bins.append({
        "label": f"ge{thresholds[-1]}",
        "desc": f">={thresholds[-1]}bp",
        "range": (thresholds[-1], float("inf")),
    })

    return bins


def classify_length(length, bins):
    """Return bin index for given length."""
    for i, bin_info in enumerate(bins):
        min_len, max_len = bin_info["range"]
        if min_len <= length < max_len:
            return i
    return len(bins) - 1  # last bin for >= max threshold


def split_bed_by_length(input_path):
    print(f"Processing {input_path} ...")

    bins = generate_bin_info(THRESHOLDS)

    if BIN_LABELS:
        if len(BIN_LABELS) != len(bins):
            raise ValueError(
                f"BIN_LABELS length ({len(BIN_LABELS)}) must equal "
                f"number of bins ({len(bins)})"
            )
        for i, label in enumerate(BIN_LABELS):
            bins[i]["label"] = label

    file_handles = []
    for bin_info in bins:
        label = bin_info["label"]
        fh = open(f"{OUT_PREFIX}_{label}.bed", "w")
        file_handles.append(fh)

    counts = [0] * len(bins)
    skipped = 0

    with open(input_path, "r") as f:
        for line in f:
            if not line.strip() or line.startswith("#") or line.startswith("track"):
                continue

            cols = line.strip().split()
            if len(cols) < 12:
                skipped += 1
                continue

            block_sizes_str = cols[10]
            try:
                block_sizes = [
                    int(x) for x in block_sizes_str.strip(",").split(",") if x
                ]
                transcript_length = sum(block_sizes)
                bin_idx = classify_length(transcript_length, bins)
                file_handles[bin_idx].write(line)
                counts[bin_idx] += 1
            except ValueError:
                skipped += 1
                continue

    for fh in file_handles:
        fh.close()

    print("Done!")
    for i, bin_info in enumerate(bins):
        print(f"{bin_info['desc']}: {counts[i]} transcripts")
    print(f"Skipped/Error: {skipped}")
    print("\nOutput files generated:")
    for bin_info in bins:
        print(f"- {OUT_PREFIX}_{bin_info['label']}.bed")


if __name__ == "__main__":
    # Prefer running in the RSeQC reference directory so outputs land next to input
    # cd /path/to/reference/RSeQC_ref/
    split_bed_by_length(INPUT_BED)
