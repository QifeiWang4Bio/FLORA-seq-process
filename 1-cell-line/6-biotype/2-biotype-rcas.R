#!/usr/bin/env Rscript
# =============================================================================
# Module : 6-biotype / RCAS gene-biotype fraction from BED intervals
# Source : 1_Cell_line/0-bin/8-biotype-rcas.R
#
# Usage:
#   Rscript 2-biotype-rcas.R <gtf_file> <bed_file1> [<bed_file2> ...]
#
# For each BED file, writes next to it under ./biotype/:
#   RCAS_<bed_basename>.csv
#   RCAS_<bed_basename>.pdf
#
# Counting rules:
#   - Annotate BED intervals against GTF via RCAS::queryGff
#   - Reads overlapping multiple biotypes → "Multi annotated"
#   - All *pseudogene* types merged → "Pseudogene"
#   - Reported categories:
#       protein_coding, lncRNA, Mt_gene, Mt_rRNA, Mt_tRNA,
#       Pseudogene, rRNA, Other, Multi annotated
#   - Percentages use total annotated reads as denominator
#
# Requires: RCAS, data.table, ggplot2, GenomeInfoDb, GenomicRanges
# =============================================================================

library(RCAS)
library(data.table)
library(ggplot2)
library(GenomeInfoDb)
library(GenomicRanges)

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
    cat("Usage: Rscript 2-biotype-rcas.R <gtf_file> <bed_file1> [<bed_file2> ...]\n")
    quit(status = 1)
}

gtf_file <- args[1]
bed_files <- args[-1]

# Categories kept in the final table / figure
MAIN_TYPES <- c(
    "protein_coding", "lncRNA", "Mt_gene", "Mt_rRNA", "Mt_tRNA",
    "Pseudogene", "rRNA"
)
FINAL_ORDER <- c(
    "protein_coding", "lncRNA", "Mt_gene", "Mt_rRNA", "Mt_tRNA",
    "Pseudogene", "rRNA", "Other", "Multi annotated"
)

# -----------------------------------------------------------------------------
# Helpers
# -----------------------------------------------------------------------------

fix_chr_prefix <- function(queryRegions, gtf_seqlevels) {
    bed_seqlevels <- seqlevels(queryRegions)
    common_chr <- intersect(bed_seqlevels, gtf_seqlevels)

    if (length(common_chr) > 0) {
        return(queryRegions)
    }

    cat("Warning: no common chromosomes; trying chr-prefix fix...\n")
    if (any(grepl("^chr", bed_seqlevels)) && !any(grepl("^chr", gtf_seqlevels))) {
        seqlevels(queryRegions) <- sub("^chr", "", bed_seqlevels)
        cat("Removed 'chr' prefix from BED\n")
    } else if (!any(grepl("^chr", bed_seqlevels)) && any(grepl("^chr", gtf_seqlevels))) {
        seqlevels(queryRegions) <- paste0("chr", bed_seqlevels)
        cat("Added 'chr' prefix to BED\n")
    } else {
        return(NULL)
    }

    if (length(intersect(seqlevels(queryRegions), gtf_seqlevels)) == 0) {
        return(NULL)
    }
    queryRegions
}

summarize_biotypes <- function(overlaps, biotype_col) {
    overlaps_select <- overlaps[, c(biotype_col, "queryIndex"), with = FALSE]
    overlaps_unique <- unique(overlaps_select)

    read_biotype_count <- overlaps_unique[, .(n_biotypes = .N), by = queryIndex]
    multi_index <- read_biotype_count[n_biotypes > 1, queryIndex]
    total_annotated <- length(unique(overlaps$queryIndex))

    cat("Multi annotated reads:", length(multi_index), "\n")
    cat("Total annotated reads:", total_annotated, "\n")

    # Keep reads with a single biotype for category counts
    overlaps_single <- overlaps_unique[!queryIndex %in% multi_index, ]
    df <- overlaps_single[, .(count = .N), by = get(biotype_col)]
    df <- as.data.frame(df)
    colnames(df) <- c("feature", "count")
    rownames(df) <- df$feature

    # Merge all *pseudogene* labels
    pseudogene_rows <- grep("pseudogene", df$feature, ignore.case = TRUE)
    if (length(pseudogene_rows) > 0) {
        pseudogene_total <- sum(df[pseudogene_rows, "count"])
        df <- df[-pseudogene_rows, , drop = FALSE]
        pseudogene_df <- data.frame(feature = "Pseudogene", count = pseudogene_total)
        rownames(pseudogene_df) <- "Pseudogene"
        df <- rbind(df, pseudogene_df)
    }

    # Multi-annotated bucket
    multi_df <- data.frame(feature = "Multi annotated", count = length(multi_index))
    rownames(multi_df) <- "Multi annotated"
    df <- rbind(df, multi_df)

    # Other = annotated reads not in main categories / multi
    main_rows <- intersect(rownames(df), MAIN_TYPES)
    counted_reads <- sum(df[c(main_rows, "Multi annotated"), "count"])
    other_count <- total_annotated - counted_reads
    other_df <- data.frame(feature = "Other", count = other_count)
    rownames(other_df) <- "Other"
    df <- rbind(df, other_df)

    existing <- intersect(FINAL_ORDER, rownames(df))
    df <- df[existing, , drop = FALSE]
    df$percent <- round(df$count / total_annotated * 100, 1)

    if ("protein_coding" %in% rownames(df)) {
        df["protein_coding", "feature"] <- "Protein_coding"
    }
    df
}

plot_biotype <- function(df, pdf_file) {
    fig <- ggplot2::ggplot(df, aes(x = reorder(feature, -percent), y = percent)) +
        geom_bar(stat = "identity", aes(fill = feature)) +
        geom_label(aes(y = percent + 0.5), label = df$percent) +
        labs(x = "Biotype", y = "Percent of annotated reads (%)") +
        theme_bw(base_size = 18) +
        theme(
            axis.text.x = element_text(size = 14, angle = 45, vjust = 0.5, hjust = 0.5),
            legend.position = "none",
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank()
        )
    ggsave(fig, file = pdf_file, width = 10, height = 6)
}

# -----------------------------------------------------------------------------
# Main
# -----------------------------------------------------------------------------

cat("Loading GTF annotation...\n")
cat("GTF file:", gtf_file, "\n")
gff <- importGtf(filePath = gtf_file, keepStandardChr = FALSE)
gtf_seqlevels <- seqlevels(gff)

for (bed_file in bed_files) {
    cat("\n==========================================\n")
    cat("Processing with RCAS:", bed_file, "\n")

    if (!file.exists(bed_file)) {
        cat("Warning: BED file not found:", bed_file, "\n")
        next
    }

    bed_dir <- dirname(bed_file)
    bed_name <- sub("\\.bed$", "", basename(bed_file))
    biotype_dir <- file.path(bed_dir, "biotype")
    if (!dir.exists(biotype_dir)) {
        dir.create(biotype_dir, recursive = TRUE)
    }
    cat("Output directory:", biotype_dir, "\n")

    queryRegions <- importBed(filePath = bed_file, keepStandardChr = FALSE)
    cat("BED chromosomes (first 5):", head(seqlevels(queryRegions), 5), "\n")
    cat("GTF chromosomes (first 5):", head(gtf_seqlevels, 5), "\n")

    queryRegions <- fix_chr_prefix(queryRegions, gtf_seqlevels)
    if (is.null(queryRegions)) {
        cat("Error: cannot reconcile chromosome names; skipping\n")
        next
    }
    cat("Found", length(intersect(seqlevels(queryRegions), gtf_seqlevels)),
        "common chromosomes\n")

    overlaps <- tryCatch({
        as.data.table(queryGff(queryRegions = queryRegions, gffData = gff))
    }, error = function(e) {
        cat("Error in queryGff:", conditionMessage(e), "\n")
        NULL
    })

    if (is.null(overlaps) || nrow(overlaps) == 0) {
        cat("Warning: no overlaps; skipping\n")
        next
    }

    biotype_col <- grep("gene_biotype", colnames(overlaps), value = TRUE)
    if (length(biotype_col) == 0) {
        cat("Error: gene_biotype column not found in GTF; skipping\n")
        next
    }

    df <- summarize_biotypes(overlaps, biotype_col)
    print(df)

    csv_file <- file.path(biotype_dir, paste0("RCAS_", bed_name, ".csv"))
    write.csv(df, file = csv_file, row.names = FALSE, quote = FALSE)
    cat("CSV saved to:", csv_file, "\n")

    pdf_file <- file.path(biotype_dir, paste0("RCAS_", bed_name, ".pdf"))
    plot_biotype(df, pdf_file)
    cat("PDF saved to:", pdf_file, "\n")
    cat("RCAS analysis completed for", bed_name, "\n")
}

cat("\n==========================================\n")
cat("All RCAS analyses completed!\n")
