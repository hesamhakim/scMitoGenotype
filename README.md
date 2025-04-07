# scMitoGenotype: Mitochondrial Genotyping from Single-Cell RNA-Seq Data

## Summary

This repository provides a collection of scripts designed to identify mitochondrial DNA (mtDNA) variants and quantify their allele frequency (AF) at the single-cell level using 10x Genomics single-cell RNA-seq (scRNA-seq) data. The workflow begins with processing raw FASTQ files using 10x Genomics Cell Ranger, followed by mtDNA variant genotyping from the resulting BAM files using one of three alternative tools: mgatk, cellSNP-lite, or SCReadCounts. Post-processing scripts are included to filter variants based on user-defined criteria (e.g., read depth, allele frequency) and generate standardized cell-by-variant tables suitable for downstream analysis.

## Workflow Overview

The analysis pipeline consists of the following major steps:

1.  **Initial scRNA-seq Processing:**
    *   Raw FASTQ files are processed using `cellranger count` (via `cellranger_RNA_count.slurm`) to generate aligned BAM files, cell barcodes, and feature matrices.

2.  **Mitochondrial Variant Genotyping (Choice of Tool):**
    *   The aligned BAM file from Cell Ranger is processed using **one** of the following tools to identify potential mtDNA variants and quantify read support per cell:
        *   **mgatk:** Utilizes `mgatk tenx` (via `mgatk_genotype.slurm`). Output is primarily a `.rds` object.
        *   **cellSNP-lite:** Utilizes `cellsnp-lite` (via `cellSNP_lite.slurm`). Output includes VCF and MTX files.
        *   **SCReadCounts:** Utilizes the `SCReadCounts` binary (via `scReadCounts.slurm`), typically run on chunks of a predefined SNV list. Output consists of multiple TSV files.

3.  **Post-Processing and Filtering:**
    *   The output from the chosen genotyping tool is processed using a corresponding R script (`mgatk_analysis.R`, `process_cellSNP_output.R`, or `process_scReadCounts_output.R`).
    *   *(SCReadCounts only)*: Output TSV files from chunks must be concatenated into a single file before running the R script.
    *   These R scripts calculate allele frequencies (if not already present), apply user-defined filters based on read depth (DP) and allele frequency (AF), and generate a final filtered table.

4.  **Final Output:**
    *   A comma-separated value (CSV) file for each sample, containing columns `variant_id`, `cell_id`, `DP`, and `AF` for variants and cells passing the specified filters.

## System Requirements

### Software Dependencies

*   **10x Genomics Cell Ranger:** (e.g., v7.x or compatible version)
*   **mgatk:** Required for the mgatk pipeline branch. (Installation via Conda recommended).
*   **cellSNP-lite:** Required for the cellSNP-lite pipeline branch. (Installation via Conda recommended).
*   **SCReadCounts:** Required for the SCReadCounts pipeline branch. (Binary executable required).
*   **SLURM Workload Manager:** Scripts are written for SLURM job submission. Adaptation may be required for other schedulers (e.g., SGE, PBS).
*   **R:** Version >= 4.0 recommended.
*   **Conda:** Recommended for managing software environments.
*   **Standard Unix/Linux command-line tools:** `awk`, `gunzip`, `cat`, `mkdir`, `bash`, etc.

### R Package Dependencies

The following R packages are required for the post-processing scripts:

*   `tidyverse`
*   `Matrix`
*   `data.table`
*   `optparse`
*   `logr`
*   `SummarizedExperiment` (from Bioconductor)
*   `vcfR`

Installation command within R:
```R
# Install CRAN packages
install.packages(c("tidyverse", "Matrix", "data.table", "optparse", "logr", "vcfR"))

# Install Bioconductor packages
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("SummarizedExperiment")
```

### Input Data

*   Raw scRNA-seq FASTQ files (typically from 10x Genomics platform).
*   A Cell Ranger compatible reference transcriptome (e.g., GRCh38).
*   *(For SCReadCounts)*: A directory containing pre-split mitochondrial SNV list files (format: `chr\tpos\tref\talt`). See `scReadCounts.slurm` configuration.

## Installation and Setup

1.  **Obtain Repository:** Clone the repository using Git:
    ```bash
    git clone https://github.com/your-username/scMitoGenotype.git
    cd scMitoGenotype
    ```

2.  **Install Software:** Ensure all software dependencies (Cell Ranger, mgatk, cellSNP-lite, SCReadCounts, R, Conda) are installed and accessible within your computing environment. Use Conda environments or module systems as appropriate.

3.  **Install R Packages:** Install the required R packages using the commands provided in the "R Package Dependencies" section.

4.  **Prepare Configuration File (`run_list.txt`):**
    *   Create a **tab-separated** text file (e.g., `my_run_list.txt`). This file defines the samples to be processed and their associated paths.
    *   Each line corresponds to one sample.
    *   The columns **must** be in the following order (or the `awk` commands within the `.slurm` scripts must be adjusted):
        1.  `ID`: A unique numeric identifier for each sample, starting from 1. This corresponds to the SLURM array task ID.
        2.  `SampleName` / `RunName`: A unique name for the sample or run, used for directory and file naming.
        3.  `FastqDir`: Full path to the directory containing FASTQ files for the sample (used by `cellranger_RNA_count.slurm`).
        4.  `TranscriptomePath`: Full path to the Cell Ranger reference transcriptome directory (used by `cellranger_RNA_count.slurm`).
        5.  `ExpectedCells`: Estimated number of cells for the sample (used by `cellranger_RNA_count.slurm`).
        6.  `ProjectDir`: The base directory path where output subdirectories (e.g., `cellranger_RNA_out`, `mgatk_out`) will be created or accessed.
    *   *Example `my_run_list.txt`:*
        ```text
        1	SampleA	/path/to/data/fastqs/SampleA	/path/to/refs/refdata-GRCh38	5000	/path/to/project/output
        2	SampleB	/path/to/data/fastqs/SampleB	/path/to/refs/refdata-GRCh38	6000	/path/to/project/output
        ```

5.  **Configure SLURM Scripts:**
    *   **Critically, edit each `.slurm` script** (`cellranger_RNA_count.slurm`, `mgatk_genotype.slurm`, `cellSNP_lite.slurm`, `scReadCounts.slurm`):
        *   Update the `RUN_LIST_FILE` variable to  point to the **absolute path** of your `run_list.txt` file.
        *   Modify the `#SBATCH` directives (e.g., `--qos`, `--account`, `--time`, `--mem`, `--cpus-per-task`) to match your cluster's requirements and resource availability.
        *   Set the `#SBATCH --array=1-N` directive, where `N` is the total number of samples (lines) in your `run_list.txt`.
        *   Verify and adjust the "Environment Setup" section in each script to correctly load necessary modules or activate Conda environments required for the specific tool being run.
        *   **For `scReadCounts.slurm`:** Update `SCREADCOUNTS_BIN_DIR` (path to the SCReadCounts executable directory) and `SNV_LIST_DIR` (path to the directory containing pre-split SNV files).
    *   Create a directory (e.g., `logs`) in the location from which you will submit jobs, to store SLURM output and error files, or modify the `#SBATCH -o` and `#SBATCH -e` paths accordingly.

## Usage Instructions

The workflow involves submitting SLURM jobs sequentially. Ensure each step completes successfully before proceeding to the next.

1.  **Run Cell Ranger:** Submit the Cell Ranger processing job for all samples defined in the array:
    ```bash
    sbatch cellranger_RNA_count.slurm
    ```
    Monitor the jobs using SLURM commands (e.g., `squeue`) and wait for completion. Check logs for errors.

2.  **Run Genotyping:** Choose **one** genotyping tool and submit the corresponding SLURM job:
    *   **Option A (mgatk):**
        ```bash
        sbatch mgatk_genotype.slurm
        ```
    *   **Option B (cellSNP-lite):**
        ```bash
        sbatch cellSNP_lite.slurm
        ```
    *   **Option C (SCReadCounts):**
        ```bash
        sbatch scReadCounts.slurm
        ```
    Monitor the jobs and wait for completion. Check logs for errors.

3.  **(SCReadCounts Only) Concatenate Output Files:** If SCReadCounts was used, combine the chunked output files (`scReadCounts_*.tsv`) into a single file for each sample before running the R script. This can be done using standard command-line tools.
    *   *Example for one sample:*
        ```bash
        # Define variables based on your run_list.txt and paths
        PROJECT_DIR="/path/to/project/output"
        RUN_NAME="SampleA"
        OUTPUT_DIR="${PROJECT_DIR}/scReadCounts_out/${RUN_NAME}/MT"

        echo "Concatenating output for ${RUN_NAME}..."
        cat ${OUTPUT_DIR}/scReadCounts_*.tsv | gzip > ${OUTPUT_DIR}/scReadCounts_filtered.tsv.gz
        echo "Concatenation complete for ${RUN_NAME}."
        ```
    Repeat this concatenation for each sample processed with SCReadCounts.

4.  **Run R Post-Processing:** Execute the R script corresponding to the genotyping tool used in Step 2. Provide appropriate command-line arguments for input paths, output file names, and filtering parameters. Run this step for each sample individually.
    *   *Example (mgatk):*
        ```bash
        Rscript mgatk_analysis.R \
            --rds_file /path/to/project/output/mgatk_out/SampleA/final/SampleA.rds \
            --output_prefix /path/to/project/output/results/SampleA_filtered_mgatk \
            --min_dp 10 \
            --min_af 0.03 \
            --max_af 0.97 \
            --log_file logs/SampleA_mgatk_analysis.log
            # Add other relevant parameters like --min_cells_conf_detected if needed
        ```
    *   *Example (cellSNP-lite):*
        ```bash
        Rscript process_cellSNP_output.R \
            --input_dir /path/to/project/output/scSNP_lite_out/SampleA \
            --output_file /path/to/project/output/results/SampleA_filtered_cellSNP.csv \
            --min_dp 10 \
            --min_ad 2 \
            --min_af 0.03 \
            --max_af 0.90 \
            --log_file logs/SampleA_process_cellSNP.log
        ```
    *   *Example (SCReadCounts):*
        ```bash
        Rscript process_scReadCounts_output.R \
            --input_file /path/to/project/output/scReadCounts_out/SampleA/MT/scReadCounts_filtered.tsv.gz \
            --output_file /path/to/project/output/results/SampleA_filtered_scReadCounts.csv \
            --min_dp 10 \
            --min_af 0.03 \
            --max_af 0.90 \
            --log_file logs/SampleA_process_scReadCounts.log
        ```
    Ensure output directories (e.g., `/path/to/project/output/results/`) exist before running the R scripts.

## Output Description

The primary output of this workflow is a set of CSV files (one per sample per tool used), typically saved in a results directory specified during the R script execution. Each file contains the filtered list of mitochondrial variants detected in single cells, with the following columns:

*   `variant_id`: A unique identifier for the mitochondrial variant (e.g., `A1624G`, combining reference allele, position, and alternate allele).
*   `cell_id`: The cell barcode identifying the single cell.
*   `DP`: Read Depth. The number of reads covering the variant position in that specific cell, supporting the allele frequency calculation (derived from `coverage` in mgatk, `DP` in cellSNP-lite, `GoodReads` in SCReadCounts).
*   `AF`: Allele Frequency. The fraction of reads supporting the alternate allele at the variant position in that specific cell.

Intermediate files generated by Cell Ranger and the genotyping tools are stored in subdirectories specified within the SLURM scripts (e.g., `cellranger_RNA_out`, `mgatk_out`, `scSNP_lite_out`, `scReadCounts_out`).

## License

This project is distributed under the terms of the MIT License. See the `LICENSE` file for details.

## Citation

If you use these scripts or the integrated tools in your research, please cite the original publications for the respective tools:

*   **mgatk:** Lareau, C.A., Ludwig, L.S., Muus, C. et al. Massively parallel single-cell mitochondrial DNA genotyping and chromatin profiling. *Nat Biotechnol* 39, 451–461 (2021). [https://doi.org/10.1038/s41587-020-00791-x](https://doi.org/10.1038/s41587-020-00791-x)
*   **cellSNP-lite:** Xu, R., Huang, Y., & Yuan, H. cellsnp-lite: an efficient tool for genotyping single cells. *Bioinformatics*, 37(23), 4569–4571 (2021). [https://doi.org/10.1093/bioinformatics/btab340](https://doi.org/10.1093/bioinformatics/btab340)
*   **SCReadCounts:** Tian, T., Zhang, J., Liu, J. et al. SCReadCounts: estimation of cell-level SNVs expression from scRNA-seq data. *BMC Genomics* 22, 670 (2021). [https://doi.org/10.1186/s12864-021-07974-8](https://doi.org/10.1186/s12864-021-07974-8)

```
