# Breast cancer and CAF analysis

This directory contains the FLORA-seq breast cancer workflow. The main data
objects are sample-level count matrices, CellMemory predictions, CAF subtype
annotations and STARsolo gene/splice-junction matrices.

## Scripts

| Files | Analysis goal |
|---|---|
| `1-process-caf.sh`, `2-coverage.sh` | Process reads and evaluate coverage. |
| `3-anno-pipeline.py`, `4-anno.py` | Prepare train/query objects and transfer coarse cell labels. |
| `5-caf-cluster.py`, `6-caf-analysis.py` | Identify high-confidence CAFs, cluster subtypes and perform differential/pathway analyses. |
| `7-sj-process.sh`, `8-sj-marvel.R` | Generate CAF splice-junction matrices and compare mCAF with tCAF using MARVEL. |

Files are numbered in run order. The final CAF subtypes are mCAF, iCAF, tCAF,
dCAF and vCAF; subtype labels originate from the unsupervised CAF clustering
workflow rather than a supervised subtype-transfer model.
