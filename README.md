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
