#!/usr/bin/env Rscript
# Usage: Rscript 4-biotype-percell-rcas.R <gtf> <bed> <output_csv>
suppressPackageStartupMessages({
    library(RCAS)
    library(data.table)
    library(GenomicRanges)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) stop("Usage: Rscript 4-biotype-percell-rcas.R <gtf> <bed> <output_csv>")
gtf_file <- args[1]; bed_file <- args[2]; out_csv <- args[3]

message("Loading GTF: ", gtf_file)
gff <- importGtf(filePath = gtf_file, keepStandardChr = FALSE)

message("Importing BED: ", bed_file)
qr  <- importBed(filePath = bed_file, keepStandardChr = FALSE)

# Auto-fix chr-prefix mismatch between BED and GTF
b <- seqlevels(qr); g <- seqlevels(gff)
if (length(intersect(b, g)) == 0) {
    if (any(grepl("^chr", b)) && !any(grepl("^chr", g))) {
        seqlevels(qr) <- sub("^chr", "", b)
    } else if (!any(grepl("^chr", b)) && any(grepl("^chr", g))) {
        seqlevels(qr) <- paste0("chr", b)
    } else {
        stop("Cannot reconcile chromosome naming between BED and GTF.")
    }
}

message("Running queryGff (this may take a few minutes)...")
ovl <- as.data.table(queryGff(queryRegions = qr, gffData = gff))

bt_col <- grep("gene_biotype", names(ovl), value = TRUE)[1]
if (is.na(bt_col)) stop("gene_biotype column not found in GTF")

# Priority-based biotype assignment per read (same scheme as bulk script)
priority <- c("protein_coding","lncRNA","snRNA","snoRNA","scaRNA","miRNA","misc_RNA","rRNA")
dt <- ovl[, .(queryIndex, biotype = get(bt_col))]
dt[, pri := match(biotype, priority)][is.na(pri), pri := 999L]
dt <- dt[dt[, .I[which.min(pri)], by = queryIndex]$V1]

# Map integer queryIndex → CB barcode (stored in mcols(qr)$name = BED col 4)
dt[, CB := mcols(qr)$name[queryIndex]]

# Collapse low-frequency biotypes
dt[!biotype %in% priority, biotype := "Other"]

fwrite(dt[, .(CB, biotype)], out_csv)
message(sprintf("Saved %s  (%d reads annotated, %d unique cells)",
                out_csv, nrow(dt), uniqueN(dt$CB)))
