# FLORA-seq raw FASTQ preprocessing

## Raw input data

Raw sequencing data were provided as paired-end FASTQ files. The cell-line
data set consisted of `scMix_3R_L01_R1.fq.gz` and `scMix_3R_L01_R2.fq.gz`.
For both the CAF and mouse hypothalamus (mHypo) data sets, two sequencing
runs were generated from the same experimental library to increase the
total number of sequencing reads, rather than representing independent
library preparations or experimental replicates. The CAF runs corresponded
to the file prefixes `2P221113028US2S2718DX_scB_P3_L01` and
`2P221114048US2S2699DX_scB_P3_L01`, whereas the mHypo runs corresponded to
`snHy1_R1/R2.fq.gz` and `snHy2_R1/R2.fq.gz`. Throughout this document these
two runs are still referred to as "batches" for processing purposes (both
are demultiplexed and reconstructed the same way and then merged), but this
term does not imply biological or technical replication.

## Reconstruction of Cell Ranger reads

Sequencer R1 contained the cDNA insert, whereas sequencer R2 encoded a 16-nt
cell barcode, an RT-primer index, and a 10-nt UMI. Samples and primer types
were identified by exact matching of the RT index in sequencer R2. The 3′ read
structure was `CB(16)–index(8/6)–UMI(10)`, whereas the random-primer read
structure was `CB(16)–UMI(10)–index(8/6)`. The 16-nt cell barcode was
reverse-complemented, its quality scores were reversed in the same order, and
the barcode was concatenated with the 10-nt UMI. The resulting CB–UMI read was
used as Cell Ranger R1, and the cDNA sequence from sequencer R1 was used as
Cell Ranger R2. Read pairs without a matching index were excluded.

## Cell-line data

The index `ATGCCTAA` identified 3′ reads, whereas `AACGTGAT` identified
random-primer reads. For the stand-alone 3′ input, CB/UMI reconstruction
(truncated to a 26-nt CB–UMI read) was followed by fastp filtering with
`-q 20 -u 30 -n 10 -w 4`. For the stand-alone random-primer input,
demultiplexed read pairs were filtered before CB/UMI reconstruction using
fastp with `-q 20 -u 30 -n 10 -w 4 -l 150`, followed by CB/UMI reconstruction;
this reconstruction only removed the RT index and retained the full
downstream sequence rather than truncating to 26 nt, consistent with the
158/142-nt `mixR` R1 observed in production data. For the combined input, 3′
and random-primer reads were identified separately and both were reconstructed
as 26-nt CB–UMI reads. The two types of read pairs were
combined and filtered with fastp using `-q 20 -u 30 -n 10 -w 4`. The filtered
files were formatted as the `mix3`, `mixR`, and `mixR3` Cell Ranger inputs.

## CAF data

The two sequencing batches were processed separately. Three 8-nt indices
(`AACGTGAT`, `AAACATCG`, and `ATGCCTAA`) assigned reads to samples 1–3, and
the index position distinguished 3′ from random-primer reads. The 3′ arm was
reconstructed as a 26-nt CB–UMI read. For the random-primer arm, the internal
RT index was removed while the sequence downstream of the index was retained,
producing a 158-nt R1 in the processed data. Corresponding sample- and
primer-specific read pairs from the two batches were then merged. The 3′ and
random-primer read pairs from each sample were subsequently combined to
generate the `sample*_3r` input.

## mHypo data

The two sequencing batches were demultiplexed and quality controlled
separately. Groups 1–12 used 8-nt indices, whereas groups 13–24 used 6-nt
indices. For the first batch, the 3′ CB–UMI read was truncated to 26 nt before
fastp; downstream sequence after the index was retained before quality control
for first-batch random-primer reads and for reconstructed reads from the
second batch. All 3′ and random-primer read pairs were filtered with fastp
using `-q 20 -u 30 -n 10 -w 12 -l 26 --cut_right`. Random-primer read pairs
were subsequently processed by two paired Trim Galore runs. The first used
`-q 20 --phred33 --stringency 3 --length 26 -e 0.1`; the second used the same
settings with the additional `--illumina` option. Cell Ranger R1 was then
truncated to 26 nt. Finally, corresponding group- and primer-specific read
pairs from the two sequencing batches were merged, and the 3′ and
random-primer read pairs were combined for each group to generate the 3r
input.
