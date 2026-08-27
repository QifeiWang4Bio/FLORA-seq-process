# FLORA-seq FASTQ preprocessing

This directory converts raw sequencer FASTQs into Cell Ranger-compatible
FASTQs for the Cell-line, CAF and mHypo data sets.

## Processing logic

1. Demultiplex read pairs by exact RT-index matching.
2. Reconstruct Cell Ranger R1 from the reverse-complemented 16-nt cell
   barcode and the 10-nt UMI encoded in sequencer R2.
3. Use sequencer R1, which contains cDNA, as Cell Ranger R2.
4. Apply data-set-specific processing and quality-control steps where used.
5. Merge corresponding primer types and sequencing batches.


## Run

Copy `reproducible/config.example.sh`, replace the input and output paths, and
run:

```bash
bash reproducible/01_cell_line.sh /path/to/config.sh
bash reproducible/02_caf.sh /path/to/config.sh
bash reproducible/03_mhypo.sh /path/to/config.sh
```

Cell Ranger-ready FASTQs are written to each data set's `cellranger_fastq/`
output directory.
