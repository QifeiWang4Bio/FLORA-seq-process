# Mouse hypothalamus analysis

This directory contains the FLORA-seq mouse hypothalamus workflow. The main
objects are the raw count matrix, C7/C25/C66 annotations, neuronal
subsets, splice-junction matrices and cell–cell interaction inputs.

## Modules

| Files | Analysis goal |
|---|---|
| `1-process-qc.R`–`3-biotype.R` | Perform sample QC, build count/HVG objects and summarize gene biotypes. |
| `4-anno-coarse.py`–`6-imprint-coarse.py` | Transfer C7/C25 labels and analyze lncRNA/imprinted-gene patterns. |
| `7-sj-process.sh`, `8-sj-marvel.R` | Compare neuronal splice-junction usage with MARVEL. |
| `9-anno-c66.py`–`11-imprint-c66.py` | Transfer C66 labels and analyze fine neuronal types. |
| `12-interaction-prep.py`, `13-interaction-c66.R` | Prepare and run C66 NeuronChat analyses. |
