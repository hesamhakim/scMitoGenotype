# scMitoGenotype: Mitochondrial Genotyping from Single-Cell RNA-Seq

<p align="center">
  <img src="https://img.shields.io/badge/Language-Bash%20%26%20R-blue" alt="Language: Bash & R">
  <img src="https://img.shields.io/badge/Pipeline-Snakemake%20%26%20SLURM-lightgrey" alt="Pipeline: SLURM">
  <img src="https://img.shields.io/badge/License-MIT-yellow.svg" alt="License: MIT">
  <img src="https://img.shields.io/badge/Status-Active-brightgreen" alt="Status: Active">
</p>

## Overview ✨

`scMitoGenotype` provides a streamlined workflow to identify mitochondrial DNA (mtDNA) variants and quantify their allele frequency (AF) at the single-cell level, starting from 10x Genomics single-cell RNA-seq (scRNA-seq) data.

After initial FASTQ processing with Cell Ranger, this repository offers scripts implementing **three distinct tools** for mtDNA variant genotyping directly from BAM files:

1.  🧬 **mgatk:** A robust tool specifically designed for mtDNA analysis in single cells.
2.  🔬 **cellSNP-lite:** An efficient tool for genotyping SNPs in single cells.
3.  📊 **SCReadCounts:** Estimates cell-level SNV expression from scRNA-seq data.

Each tool's pipeline includes SLURM-optimized processing scripts and corresponding R scripts for parsing, filtering (based on depth, AF), and generating standardized cell-by-variant tables.

## Workflow Diagram 🗺️

```mermaid
graph TD
    subgraph Input
        A[FASTQ Files]
        B(Reference Transcriptome)
    end

    subgraph Step 1: Cell Ranger Processing
        C{cellranger_RNA_count.slurm}
        D[Aligned BAM File]
        E[Cell Barcodes]
    end

    subgraph Step 2: Mitochondrial Genotyping (Choose One)
        direction LR
        subgraph Option A: mgatk
            F{mgatk_genotype.slurm} --> G((mgatk .rds Output))
        end
        subgraph Option B: cellSNP-lite
            H{cellSNP_lite.slurm} --> I((VCF & MTX Output))
        end
        subgraph Option C: SCReadCounts
            J{scReadCounts.slurm} --> K((Chunked TSV Output))
        end
    end

    subgraph Step 3: Post-Processing & Filtering
        direction LR
        subgraph Option A Post
             L(mgatk_analysis.R)
        end
        subgraph Option B Post
             M(process_cellSNP_output.R)
        end
        subgraph Option C Post
             N(Concatenate TSV Chunks) --> O(process_scReadCounts_output.R)
        end
    end

    subgraph Output
        P[Filtered Cell x Variant Table (.csv)]
    end

    Input --> C;
    C --> D;
    C --> E;
    D --> F; E --> F;
    D --> H; E --> H;
    D --> J; E --> J; B --> J; # SCReadCounts also needs SNV list implicitly

    G --> L;
    I --> M;
    K --> N;
    O --> P;
    L --> P;
    M --> P;

    style C fill:#f9f,stroke:#333,stroke-width:2px
    style F fill:#ccf,stroke:#333,stroke-width:2px
    style H fill:#cfc,stroke:#333,stroke-width:2px
    style J fill:#fec,stroke:#333,stroke-width:2px
    style L fill:#ccf,stroke:#333,stroke-width:1px,stroke-dasharray: 5 5
    style M fill:#cfc,stroke:#333,stroke-width:1px,stroke-dasharray: 5 5
    style N fill:#fec,stroke:#333,stroke-width:1px,stroke-dasharray: 5 5
    style O fill:#fec,stroke:#333,stroke-width:1px,stroke-dasharray: 5 5
    style P fill:#bbf,stroke:#333,stroke-width:2px

```

## Key Features 🚀

*   **Flexible Genotyping:** Choose from three popular tools (mgatk, cellSNP-lite, SCReadCounts).
*   **SLURM Integration:** Scripts are optimized for SLURM cluster environments using array jobs for easy parallelization across samples.
*   **Configurable Filtering:** R scripts allow easy adjustment of depth (DP) and allele frequency (AF) thresholds.
*   **Standardized Output:** Generates consistent CSV tables (`variant_id`, `cell_id`, `DP`, `AF`) regardless of the tool used.
*   **Reproducibility:** Provides a clear, script-based workflow from FASTQ to variant tables.

## Repository Contents 📁

*   `cellranger_RNA_count.slurm`: Runs `cellranger count`.
*   `mgatk_genotype.slurm`: Runs `mgatk tenx`.
*   `mgatk_analysis.R`: Processes `mgatk` output.
*   `cellSNP_lite.slurm`: Runs `cellsnp-lite`.
*   `process_cellSNP_output.R`: Processes `cellsnp-lite` output.
*   `scReadCounts.slurm`: Runs `SCReadCounts` (in chunks).
*   `process_scReadCounts_output.R`: Processes concatenated `SCReadCounts` output.
*   `run_list.txt`: **Example** configuration file. **Create your own!**
*   `README.md`: This file.
*   `LICENSE`: Project license file (e.g., MIT).

## Prerequisites 🛠️

### Software

*   **10x Genomics Cell Ranger:** (e.g., v7.x)
*   **mgatk:** (Install via Conda recommended)
*   **cellSNP-lite:** (Install via Conda recommended)
*   **SCReadCounts:** (Binary or compiled)
*   **SLURM Workload Manager**
*   **R:** (>= 4.0 recommended)
*   **Conda:** (Recommended for environment management)
*   **Standard Unix tools:** `awk`, `gunzip`, `cat`, `mkdir`, etc.

### R Packages

```R
# Run in R console
install.packages(c("tidyverse", "Matrix", "data.table", "optparse", "logr", "vcfR"))

if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("SummarizedExperiment")
```

### Input Data

*   Raw scRNA-seq FASTQ files (10x Genomics).
*   Cell Ranger compatible reference transcriptome.
*   (For SCReadCounts) Pre-split mitochondrial SNV list files.

---

## Quick Start Instructions 🏁

1.  **Setup ⚙️:**
    *   Clone repo: `git clone ...`
    *   Install prerequisites (Software & R packages).
    *   Create **your** `run_list.txt` (tab-separated: `ID`, `SampleName`, `FastqDir`, `TranscriptomePath`, `ExpectedCells`, `ProjectDir`).
    *   Edit **all** `.slurm` scripts: Update `RUN_LIST_FILE` path, `#SBATCH` directives (qos, account, resources, array range), and environment setup (modules/conda).
    *   Create `logs` directory: `mkdir logs`.

2.  **Usage ▶️:**
    *   Run Cell Ranger: `sbatch cellranger_RNA_count.slurm`. Wait.
    *   Choose **one** genotyping tool (mgatk, cellSNP-lite, or SCReadCounts).
    *   Run chosen tool's SLURM script: `sbatch mgatk_genotype.slurm` **OR** `sbatch cellSNP_lite.slurm` **OR** `sbatch scReadCounts.slurm`. Wait.
    *   **(SCReadCounts Only) Concatenate Output 🧩:** Combine chunked `.tsv` files: `cat path/to/scReadCounts_*.tsv | gzip > path/to/scReadCounts_filtered.tsv.gz`.
    *   Run corresponding R script (e.g., `Rscript mgatk_analysis.R --rds_file ... --output_prefix ...`) with correct paths and filtering parameters.

---

## Detailed Setup 🔧

1.  **Clone Repository:**
    ```bash
    git clone https://github.com/your-username/scMitoGenotype.git
    cd scMitoGenotype
    ```

2.  **Install Dependencies:**
    *   Ensure Cell Ranger, mgatk, cellSNP-lite, and SCReadCounts are installed and accessible (e.g., via `module load` or Conda). See Prerequisites.
    *   Install required R packages (see Prerequisites).
    *   *Example Conda setup:*
        ```bash
        conda create -n mgatk_env mgatk picard R -c bioconda -c conda-forge
        conda create -n cellsnp_env cellsnp-lite R -c bioconda -c conda-forge
        # Ensure SCReadCounts binary is accessible
        ```

3.  **Configure `run_list.txt`:**
    *   Create a **tab-separated** file (e.g., `my_run_list.txt`). Columns **must** match the order expected by the `awk` commands in the `.slurm` scripts (or modify the scripts):
        1.  `ID`: Numeric ID (1, 2, 3...) matching SLURM array index.
        2.  `SampleName` / `RunName`: Sample identifier.
        3.  `FastqDir`: Path to FASTQ directory (for `cellranger_RNA_count.slurm`).
        4.  `TranscriptomePath`: Path to Cell Ranger reference (for `cellranger_RNA_count.slurm`).
        5.  `ExpectedCells`: Estimated cells (for `cellranger_RNA_count.slurm`).
        6.  `ProjectDir`: Base path for input/output subdirectories.
    *   *Example `my_run_list.txt`:*
        ```text
        1	SampleA	/path/to/fastqs/SampleA	/path/to/refdata-GRCh38	5000	/project/base/path
        2	SampleB	/path/to/fastqs/SampleB	/path/to/refdata-GRCh38	6000	/project/base/path
        ```

4.  **Configure SLURM Scripts:**
    *   **Crucially, edit each `.slurm` script:**
        *   Set `RUN_LIST_FILE` to your `run_list.txt` path.
        *   Adjust `#SBATCH` directives (qos, account, time, mem, cpus, array range `1-N`).
        *   Verify/update the "Environment Setup" section (modules/conda activation).
        *   **`scReadCounts.slurm`:** Update `SCREADCOUNTS_BIN_DIR` and `SNV_LIST_DIR`.
    *   Create `logs` directory: `mkdir logs` (or update log paths in scripts).

## Detailed Usage 🚀

Execute sequentially for each sample array job.

1.  **Run Cell Ranger:**
    ```bash
    sbatch cellranger_RNA_count.slurm
    ```
    Monitor jobs and wait for successful completion.

2.  **Run Genotyping (Choose ONE tool):**
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
    Monitor jobs and wait for successful completion.

3.  **(SCReadCounts Only) Concatenate Output:**
    *   Combine the output chunks for each sample.
    *   *Example:*
        ```bash
        PROJECT_DIR="/project/base/path" # From run_list.txt
        RUN_NAME="SampleA"             # From run_list.txt
        OUTPUT_DIR="${PROJECT_DIR}/scReadCounts_out/${RUN_NAME}/MT"
        echo "Concatenating SCReadCounts output for ${RUN_NAME}..."
        cat ${OUTPUT_DIR}/scReadCounts_*.tsv | gzip > ${OUTPUT_DIR}/scReadCounts_filtered.tsv.gz
        echo "Done."
        ```

4.  **Run R Processing Script:**
    *   Execute the R script corresponding to the chosen tool, providing appropriate arguments.
    *   *Example (mgatk):*
        ```bash
        Rscript mgatk_analysis.R \
            --rds_file /project/base/path/mgatk_out/SampleA/final/SampleA.rds \
            --output_prefix results/SampleA_filtered_mgatk \
            --min_dp 10 --min_af 0.03 --max_af 0.97 \
            --log_file logs/SampleA_mgatk_analysis.log
        ```
    *   *Example (cellSNP-lite):*
        ```bash
        Rscript process_cellSNP_output.R \
            --input_dir /project/base/path/scSNP_lite_out/SampleA \
            --output_file results/SampleA_filtered_cellSNP.csv \
            --min_dp 10 --min_ad 2 --min_af 0.03 --max_af 0.90 \
            --log_file logs/SampleA_process_cellSNP.log
        ```
    *   *Example (SCReadCounts):*
        ```bash
        Rscript process_scReadCounts_output.R \
            --input_file /project/base/path/scReadCounts_out/SampleA/MT/scReadCounts_filtered.tsv.gz \
            --output_file results/SampleA_filtered_scReadCounts.csv \
            --min_dp 10 --min_af 0.03 --max_af 0.90 \
            --log_file logs/SampleA_process_scReadCounts.log
        ```
    *   Ensure output directories (e.g., `results/`) exist. Adjust parameters as needed.

## Output Files 📊

*   **Intermediate:** Cell Ranger outputs, tool-specific outputs (`.rds`, VCF/MTX, `.tsv` chunks).
*   **Final:** A CSV file per sample/tool (e.g., `results/SampleA_filtered_mgatk.csv`) with columns:
    *   `variant_id`: Variant identifier (e.g., `A1624G`).
    *   `cell_id`: Cell barcode.
    *   `DP`: Read depth for the variant in the cell.
    *   `AF`: Allele frequency for the variant in the cell.

## Contributing 🤝

Contributions are welcome! Please feel free to submit issues or pull requests for improvements, bug fixes, or new features.

## License 📜

This project is licensed under the [MIT License](LICENSE).

## Citation 🙏

If you use these scripts or the integrated tools in your research, please cite the original publications:

*   **mgatk:** Lareau, C.A., et al. *Nat Biotechnol* 39, 451–461 (2021). [doi:10.1038/s41587-020-00791-x](https://doi.org/10.1038/s41587-020-00791-x)
*   **cellSNP-lite:** Xu, R., et al. *Bioinformatics*, 37(23), 4569–4571 (2021). [doi:10.1093/bioinformatics/btab340](https://doi.org/10.1093/bioinformatics/btab340)
*   **SCReadCounts:** Tian, T., et al. *BMC Genomics* 22, 670 (2021). [doi:10.1186/s12864-021-07974-8](https://doi.org/10.1186/s12864-021-07974-8)

```
