# Cell-line benchmark

This directory contains the FLORA-seq species-mixing and cross-platform
benchmark analyses. The main data objects are human/mouse Cell Ranger count
matrices, HEK293T/NIH3T3 annotations, depth-matched BAM files and public
SMART-seq3, VASA-seq and 10x datasets.

## Modules

| Directory | Analysis goal |
|---|---|
| `1-mapping/` | Map FLORA-seq and comparison datasets. |
| `2-genebody/` | Compare gene-body coverage across methods. |
| `3-reads-per-cell/` | Measure genes and molecules detected with increasing depth. |
| `4-sj/` | Compare splice-junction detection and validate junctions with long reads. |
| `5-matrix/` | Perform species-mixing QC and export annotated count matrices. |
| `6-biotype/` | Quantify read and count fractions by gene biotype. |

Files within each module are numbered in approximate run order. Configure
external inputs with the repository-level path settings before running.
