# Single-Sample GSEA Analysis Pipeline (SSGSEA)

## Overview
This repository contains knowhow i performed single-sample Gene Set Enrichment Analysis (ssGSEA) on epithelial cell RNA sequencing data using MSigDB hallmark gene sets and custom cancer-related gene signatures.

## Dependencies
Required R packages:
- GSVA (for ssGSEA analysis)
- Seurat (for single-cell data handling)
- ggplot2 (for visualization)
- GSEABase (for gene set handling)

Install packages using:
```R
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(c("GSVA", "GSEABase"))
install.packages(c("Seurat", "ggplot2"))
```

## Input Data
The pipeline requires:
1. A Seurat object containing epithelial cell RNA-seq data
2. MSigDB hallmark gene sets (H.gsea.rds)
3. Custom cancer-related gene signatures (provided in the code)

## File Structure
```
.
├── pdac_pancreas_epithelial.rds    # Input Seurat object
├── H.gsea.rds                      # MSigDB hallmark gene sets
└── ssGSEA_scores_epi_msigdb.rds    # Output ssGSEA scores
```

## Workflow
1. **Data Loading**
   - Loads pre-processed Seurat object
   - Imports MSigDB hallmark gene sets

2. **Gene Set Preparation**
   - Converts GeneSetCollection to list format
   - Extracts gene IDs and set names

3. **Expression Data Extraction**
   - Retrieves RNA count matrix from Seurat object

4. **ssGSEA Analysis**
   - Performs ssGSEA using GSVA package
   - Calculates enrichment scores for each gene set

5. **Results Storage**
   - Saves ssGSEA scores to RDS file

## Usage
```R
# Load required libraries
library(GSVA)
library(Seurat)
library(ggplot2)
library(GSEABase)

# Run the analysis
source("ssgsea_analysis.R")
```

## Output
The script generates "ssGSEA_scores_epi_msigdb.rds" containing enrichment scores for each gene set across all samples.


## Notes to myself
- Ensure sufficient computational resources for large datasets : Cehcked
- Default parameters for ssGSEA are used : Checked


## METHODS 

# Single-Sample Gene Set Enrichment Analysis of Epithelial Cell Transcriptomes

## Data Preprocessing and Quality Control
Single-cell RNA sequencing data from pancreatic epithelial cells were processed using the Seurat pipeline . Raw count matrices were normalized according to standard Seurat workflows. Quality control metrics were applied to remove low-quality cells and potential doublets.

Gene Set Collection
Two complementary gene set collections were utilized for the enrichment analysis:

1. MSigDB Hallmark Gene Sets: Curated hallmark gene sets were obtained from the Molecular Signatures Database (MSigDB), representing well-defined biological states and processes.

Single-Sample GSEA Implementation
Single-sample Gene Set Enrichment Analysis (ssGSEA) was performed using the GSVA package in R. The analysis was conducted using the following approach:

1. Expression data were extracted from the normalized RNA count matrix maintained in the Seurat object.

2. Gene sets were preprocessed using the GSEABase package to ensure proper formatting and gene identifier compatibility.

3. ssGSEA was implemented using the default GSVA parameters:
   - Random sample permutations for significance estimation
   - Gaussian-distributed scores
   - No minimum gene set size filter
   - No maximum gene set size filter

The ssGSEA algorithm calculated separate enrichment scores for each sample-gene set pair, enabling cell-specific pathway activity assessment.

Statistical Analysis
Enrichment scores were computed as a function of the empirical cumulative distribution functions of gene expression ranks inside and outside each gene set. This approach allows for the quantification of gene set enrichment at the single-cell level while accounting for differences in gene set sizes and expression distribution patterns.


Software Specifications
All analyses were performed in R version (?). Key package versions:
- GSVA: 
- Seurat:
- GSEABase:
- ggplot2:

References
1. Hänzelmann S, Castelo R, Guinney J. GSVA: gene set variation analysis for microarray and RNA-seq data. BMC Bioinformatics. 2013
2. Stuart T, Butler A, Hoffman P, et al. Comprehensive Integration of Single-Cell Data. Cell. 2019
3. Liberzon A, Birger C, Thorvaldsdóttir H, et al. The Molecular Signatures Database Hallmark Gene Set Collection. Cell Systems. 2015

