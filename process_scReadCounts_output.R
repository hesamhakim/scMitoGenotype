#!/usr/bin/env Rscript

# --- Purpose ---
# Processes the concatenated output file from SCReadCounts for a single sample.
# Assumes chunked outputs from the SLURM script have been combined into one file.
# Filters variants based on read counts (GoodReads/DP) and allele frequency (AF).
# Outputs a filtered table of cell-variant information.

# --- Libraries ---
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(data.table)) # Using fread for potentially large files
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(logr))

# --- Argument Parsing ---
option_list <- list(
  make_option(c("-i", "--input_file"), type = "character", default = NULL,
              help = "Path to the (concatenated) SCReadCounts output TSV file (can be gzipped)", metavar = "character"),
  make_option(c("-o", "--output_file"), type = "character", default = "filtered_variants_scReadCounts.csv",
              help = "Path for the output CSV file [default= %default]", metavar = "character"),
  make_option(c("--min_dp"), type = "integer", default = 10,
              help = "Minimum 'GoodReads' (renamed to DP) per cell-variant [default= %default]", metavar = "integer"),
  make_option(c("--min_af"), type = "double", default = 0.03,
              help = "Minimum allele frequency (AF/VAF) per cell-variant [default= %default]", metavar = "double"),
  make_option(c("--max_af"), type = "double", default = 0.90,
              help = "Maximum allele frequency (AF/VAF) per cell-variant [default= %default]", metavar = "double"),
  make_option(c("-l", "--log_file"), type = "character", default = NULL,
              help = "Optional file path for logging output.", metavar = "character")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# --- Setup Logging ---
if (!is.null(opt$log_file)) {
    log_open(opt$log_file, logdir = FALSE, show_notes = FALSE)
} else {
    log_open(stdout()) # Log to console if no file specified
}

# --- Input Validation ---
if (is.null(opt$input_file)) {
  log_print("Error: Input file (--input_file) is required.", hide_notes = TRUE)
  print_help(opt_parser)
  stop("Input file not specified.", call. = FALSE)
}
if (!file.exists(opt$input_file)) {
  log_print(paste("Error: Input file not found:", opt$input_file), hide_notes = TRUE)
  stop("Input file not found.", call. = FALSE)
}

log_print("--- Starting SCReadCounts Output Processing ---", hide_notes = TRUE)
log_print(paste("Input File:", opt$input_file), hide_notes = TRUE)
log_print(paste("Output File:", opt$output_file), hide_notes = TRUE)
log_print("Parameters:", hide_notes = TRUE)
log_print(paste("  min_dp (GoodReads):", opt$min_dp), hide_notes = TRUE)
log_print(paste("  min_af:", opt$min_af), hide_notes = TRUE)
log_print(paste("  max_af:", opt$max_af), hide_notes = TRUE)

# --- Main Workflow ---

# 1. Read SCReadCounts output file
log_print("Reading SCReadCounts output file...", hide_notes = TRUE)
variant_tbl <- tryCatch({
    fread(opt$input_file, header = FALSE, sep = "\t",
          col.names = c("CHROM", "POS", "REF", "ALT", "cell_id", "SNVCountForward", "SNVCountReverse",
                        "RefCountForward", "RefCountReverse", "SNVCount", "RefCount", "GoodReads", "PercentBadRead", "AF"))
}, error = function(e) {
    log_print(paste("Error reading input file:", e$message), hide_notes = TRUE)
    stop("Failed to read input file.", call. = FALSE)
})

log_print(paste("Read", nrow(variant_tbl), "entries from input file."), hide_notes = TRUE)

# 2. Create variant ID
log_print("Creating variant IDs...", hide_notes = TRUE)
# Ensure correct pasting order: REF, POS, ALT
variant_tbl[, variant_id := paste0(REF, POS, ALT)]

# 3. Apply Filters and select columns
log_print("Applying filters (GoodReads/DP, AF range)...", hide_notes = TRUE)
variant_tbl_filtered <- variant_tbl %>%
  filter(GoodReads >= opt$min_dp & AF >= opt$min_af & AF <= opt$max_af) %>%
  select(variant_id, cell_id, GoodReads, AF) %>%
  rename(DP = GoodReads) # Rename GoodReads to DP for consistency

log_print(paste("Final filtered table contains", nrow(variant_tbl_filtered), "cell-variant entries."), hide_notes = TRUE)

# 4. Save the final filtered table
log_print(paste("Saving filtered cell-variant table to:", opt$output_file), hide_notes = TRUE)
tryCatch({
    fwrite(variant_tbl_filtered, file = opt$output_file, row.names = FALSE, quote = FALSE, sep = ",")
}, error = function(e) {
    log_print(paste("Error writing output file:", e$message), hide_notes = TRUE)
    stop("Failed to write output file.", call. = FALSE)
})

log_print("--- Processing Complete ---", hide_notes = TRUE)
log_close()