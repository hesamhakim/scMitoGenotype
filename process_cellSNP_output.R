#!/usr/bin/env Rscript

# --- Purpose ---
# Processes the output of cellsnp-lite (.vcf, .mtx files) for a single sample.
# Calculates allele frequency (AF) per cell per variant.
# Filters variants based on depth (DP), alternate allele depth (AD), and AF.
# Outputs a filtered table of cell-variant information.

# --- Libraries ---
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(Matrix))
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(logr))
suppressPackageStartupMessages(library(vcfR)) # Using vcfR for robust VCF parsing

# --- Argument Parsing ---
option_list <- list(
  make_option(c("-i", "--input_dir"), type = "character", default = NULL,
              help = "Path to the cellsnp-lite output directory containing .vcf and .mtx files", metavar = "character"),
  make_option(c("-o", "--output_file"), type = "character", default = "filtered_variants_cellSNP.csv",
              help = "Path for the output CSV file [default= %default]", metavar = "character"),
  make_option(c("--min_dp"), type = "integer", default = 10,
              help = "Minimum read depth (DP) per cell-variant [default= %default]", metavar = "integer"),
  make_option(c("--min_ad"), type = "integer", default = 2,
              help = "Minimum alternate allele depth (AD) per cell-variant [default= %default]", metavar = "integer"),
  make_option(c("--min_af"), type = "double", default = 0.03,
              help = "Minimum allele frequency (AF) per cell-variant [default= %default]", metavar = "double"),
  make_option(c("--max_af"), type = "double", default = 0.90,
              help = "Maximum allele frequency (AF) per cell-variant [default= %default]", metavar = "double"),
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
if (is.null(opt$input_dir)) {
  log_print("Error: Input directory (--input_dir) is required.", hide_notes = TRUE)
  print_help(opt_parser)
  stop("Input directory not specified.", call. = FALSE)
}
if (!dir.exists(opt$input_dir)) {
  log_print(paste("Error: Input directory not found:", opt$input_dir), hide_notes = TRUE)
  stop("Input directory not found.", call. = FALSE)
}

# Define expected input file names
vcf_file <- file.path(opt$input_dir, "cellSNP.base.vcf.gz") # Assume gzipped from --gzip flag
dp_mtx_file <- file.path(opt$input_dir, "cellSNP.tag.DP.mtx")
ad_mtx_file <- file.path(opt$input_dir, "cellSNP.tag.AD.mtx")
# oth_mtx_file <- file.path(opt$input_dir, "cellSNP.tag.OTH.mtx") # OTH not needed for standard VAF calculation
samples_file <- file.path(opt$input_dir, "cellSNP.samples.tsv")

# Check existence of required files
required_files <- c(vcf_file, dp_mtx_file, ad_mtx_file, samples_file)
files_exist <- file.exists(required_files)
if (!all(files_exist)) {
    missing_files <- required_files[!files_exist]
    log_print("Error: Missing required input files in input directory:", hide_notes = TRUE)
    log_print(paste("  -", missing_files, collapse = "\n"), hide_notes = TRUE)
    stop("Missing input files.", call. = FALSE)
}

log_print("--- Starting cellsnp-lite Output Processing ---", hide_notes = TRUE)
log_print(paste("Input Directory:", opt$input_dir), hide_notes = TRUE)
log_print(paste("Output File:", opt$output_file), hide_notes = TRUE)
log_print("Parameters:", hide_notes = TRUE)
log_print(paste("  min_dp:", opt$min_dp), hide_notes = TRUE)
log_print(paste("  min_ad:", opt$min_ad), hide_notes = TRUE)
log_print(paste("  min_af:", opt$min_af), hide_notes = TRUE)
log_print(paste("  max_af:", opt$max_af), hide_notes = TRUE)

# --- Main Workflow ---

# 1. Read VCF to get variant information
log_print("Reading VCF file...", hide_notes = TRUE)
vcf <- tryCatch({
    read.vcfR(vcf_file, verbose = FALSE)
}, error = function(e) {
    log_print(paste("Error reading VCF file:", e$message), hide_notes = TRUE)
    stop("Failed to read VCF.", call. = FALSE)
})

variant_list <- as.data.frame(getFIX(vcf, getINFO = TRUE)) %>%
  rownames_to_column(var = "vcf_idx") %>% # Keep original row order index if needed
  select(CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO) %>%
  mutate(POS = as.integer(POS)) %>%
  mutate(variant_id = paste0(REF, POS, ALT)) %>%
  select(POS, REF, ALT, variant_id, INFO) # Keep only necessary columns

log_print(paste("Read", nrow(variant_list), "variants from VCF."), hide_notes = TRUE)

# 2. Read cell barcodes
log_print("Reading sample/cell list...", hide_notes = TRUE)
cell_ids <- read.table(samples_file, header = FALSE, stringsAsFactors = FALSE)$V1
log_print(paste("Read", length(cell_ids), "cell IDs."), hide_notes = TRUE)

# 3. Read MTX files
log_print("Reading DP matrix...", hide_notes = TRUE)
variant_tbl_DP <- tryCatch({
    Matrix::readMM(dp_mtx_file)
}, error = function(e) {
    log_print(paste("Error reading DP matrix:", e$message), hide_notes = TRUE)
    stop("Failed to read DP matrix.", call. = FALSE)
})

log_print("Reading AD matrix...", hide_notes = TRUE)
variant_tbl_AD <- tryCatch({
    Matrix::readMM(ad_mtx_file)
}, error = function(e) {
    log_print(paste("Error reading AD matrix:", e$message), hide_notes = TRUE)
    stop("Failed to read AD matrix.", call. = FALSE)
})

# Check dimensions
if (nrow(variant_tbl_DP) != nrow(variant_list) || ncol(variant_tbl_DP) != length(cell_ids)) {
    log_print("Error: DP matrix dimensions do not match VCF variants or cell IDs.", hide_notes = TRUE)
    stop("DP matrix dimension mismatch.")
}
if (nrow(variant_tbl_AD) != nrow(variant_list) || ncol(variant_tbl_AD) != length(cell_ids)) {
    log_print("Error: AD matrix dimensions do not match VCF variants or cell IDs.", hide_notes = TRUE)
    stop("AD matrix dimension mismatch.")
}

# Assign dimnames (important for converting to long format correctly)
rownames(variant_tbl_DP) <- variant_list$variant_id
colnames(variant_tbl_DP) <- cell_ids
rownames(variant_tbl_AD) <- variant_list$variant_id
colnames(variant_tbl_AD) <- cell_ids

# 4. Convert sparse matrices to long data tables
log_print("Converting matrices to long format...", hide_notes = TRUE)
summary_DP <- summary(variant_tbl_DP)
summary_AD <- summary(variant_tbl_AD)

dt_DP <- data.table(
    variant_id = rownames(variant_tbl_DP)[summary_DP$i],
    cell_id = colnames(variant_tbl_DP)[summary_DP$j],
    DP = summary_DP$x
)
dt_AD <- data.table(
    variant_id = rownames(variant_tbl_AD)[summary_AD$i],
    cell_id = colnames(variant_tbl_AD)[summary_AD$j],
    AD = summary_AD$x
)

# 5. Merge DP and AD tables and calculate AF
log_print("Merging DP and AD tables and calculating AF...", hide_notes = TRUE)
# Perform a full join to capture all cell-variant pairs with non-zero DP or AD
# Then filter out those with DP=0, as AF is undefined.
variant_tbl <- merge(dt_DP, dt_AD, by = c("variant_id", "cell_id"), all = TRUE)

# Fill NAs introduced by the merge (e.g., if AD=0, it might be missing from dt_AD)
variant_tbl[is.na(DP), DP := 0]
variant_tbl[is.na(AD), AD := 0]

# Calculate Allele Frequency (AF = AD / DP)
# Add a small epsilon to DP to prevent division by zero, although DP should ideally be >0 after filtering
variant_tbl[, AF := ifelse(DP > 0, AD / DP, 0)]

log_print(paste("Generated table with", nrow(variant_tbl), "cell-variant entries (DP > 0 or AD > 0)."), hide_notes = TRUE)

# 6. Apply Filters
log_print("Applying filters (DP, AD, AF range)...", hide_notes = TRUE)
variant_tbl_filtered <- variant_tbl %>%
  filter(DP >= opt$min_dp & AD >= opt$min_ad & AF >= opt$min_af & AF <= opt$max_af) %>%
  select(variant_id, cell_id, DP, AF) # Keep only essential columns

log_print(paste("Final filtered table contains", nrow(variant_tbl_filtered), "cell-variant entries."), hide_notes = TRUE)

# 7. Save the final filtered table
log_print(paste("Saving filtered cell-variant table to:", opt$output_file), hide_notes = TRUE)
tryCatch({
    fwrite(variant_tbl_filtered, file = opt$output_file, row.names = FALSE, quote = FALSE, sep = ",")
}, error = function(e) {
    log_print(paste("Error writing output file:", e$message), hide_notes = TRUE)
    stop("Failed to write output file.", call. = FALSE)
})

log_print("--- Processing Complete ---", hide_notes = TRUE)
log_close()