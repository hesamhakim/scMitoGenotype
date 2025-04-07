# scMitoGenotype: Mitochondrial Genotyping from Single-Cell RNA-Seq

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT) <!-- Replace MIT with your actual license if different -->

## Overview

`scMitoGenotype` provides a collection of scripts to identify mitochondrial DNA (mtDNA) variants and quantify their allele frequency (AF) at the single-cell level from 10x Genomics single-cell RNA-seq (scRNA-seq) data.

The workflow starts with standard scRNA-seq FASTQ processing using Cell Ranger. Subsequently, it offers three alternative tools for mtDNA variant genotyping directly from the Cell Ranger output BAM files:

1.  **mgatk:** (Lareau et al., 2021, *Nat Biotechnol*) - A robust tool specifically designed for mtDNA analysis in single cells.
2.  **cellSNP-lite:** (Xu et al., 2021, *Bioinformatics*) - An efficient tool for genotyping SNPs in single cells.
3.  **SCReadCounts:** (Tian et al., 2021, *BMC Genomics*) - Estimates cell-level SNV expression from scRNA-seq.

Each tool has its own processing script (optimized for SLURM environments) and a corresponding R script to parse the output, filter variants based on user-defined criteria (depth, allele frequency), and generate a standardized cell-by-variant table.

## Workflow

The general workflow is as follows:

1.  **FASTQ Processing:**
    *   `cellranger_RNA_count.slurm`: Run `cellranger count` on raw FASTQ files to generate aligned BAM files and cell barcodes.
2.  **Mitochondrial Genotyping (Choose ONE per sample):**
    *   **Option A (mgatk):**
        *   `mgatk_genotype.slurm`: Run `mgatk tenx` on the Cell Ranger BAM file.
        *   `mgatk_analysis.R`: Process the resulting `.rds` file to call variants and generate a filtered table.
    *   **Option B (cellSNP-lite):**
        *   `cellSNP_lite.slurm`: Run `cellsnp-lite` on the Cell Ranger BAM file.
        *   `process_cellSNP_output.R`: Process the resulting VCF and MTX files to generate a filtered table.
    *   **Option C (SCReadCounts):**
        *   `scReadCounts.slurm`: Run `SCReadCounts` (in chunks) on the Cell Ranger BAM file using pre-defined SNV lists.
        *   *(Manual Step)*: Concatenate the chunked output TSV files.
        *   `process_scReadCounts_output.R`: Process the concatenated TSV file to generate a filtered table.
3.  **Output:**
    *   A final CSV table for each sample and chosen tool, containing `variant_id`, `cell_id`, `DP` (Depth), and `AF` (Allele Frequency) for variants passing filters.

## Repository Contents

*   `cellranger_RNA_count.slurm`: SLURM script for running `cellranger count`.
*   `mgatk_genotype.slurm`: SLURM script for running `mgatk tenx`.
*   `mgatk_analysis.R`: R script to process `mgatk` output (`.rds`) and generate a filtered variant table.
*   `cellSNP_lite.slurm`: SLURM script for running `cellsnp-lite`.
*   `process_cellSNP_output.R`: R script to process `cellsnp-lite` output (`.vcf`, `.mtx`) and generate a filtered variant table.
*   `scReadCounts.slurm`: SLURM script for running `SCReadCounts` (in chunks).
*   `process_scReadCounts_output.R`: R script to process concatenated `SCReadCounts` output (`.tsv`) and generate a filtered variant table.
*   `run_list.txt`: **Example** configuration file defining samples, paths, and parameters for SLURM jobs. **Users must create and configure their own.**
*   `README.md`: This file.
*   `LICENSE`: Project license file (e.g., MIT).

## Prerequisites

### Software

*   **10x Genomics Cell Ranger:** For initial FASTQ processing (`cellranger count`). Version compatibility may vary; tested with vX.Y.Z.
*   **mgatk:** Required for the `mgatk` pipeline branch. Installation typically via Conda.
*   **cellSNP-lite:** Required for the `cellSNP-lite` pipeline branch. Installation typically via Conda.
*   **SCReadCounts:** Required for the `SCReadCounts` pipeline branch. Download binary or compile from source.
*   **SLURM Workload Manager:** The `.slurm` scripts are designed for clusters using SLURM. Modifications may be needed for other job schedulers.
*   **R:** Version >= 4.0 recommended.
*   **Conda:** Recommended for managing Python/tool dependencies.
*   **Standard Unix tools:** `awk`, `gunzip`, `cat`, `mkdir`, etc.

### R Packages

Install the following R packages:

```R
install.packages(c("tidyverse", "Matrix", "data.table", "optparse", "logr", "SummarizedExperiment", "vcfR"))
# Note: SummarizedExperiment is part of Bioconductor
# if (!require("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")
# BiocManager::install("SummarizedExperiment")
