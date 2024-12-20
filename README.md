# NYU Singlecell analysis
Single cell sample acquisition, processing, and analysis.

# Single-Cell RNA Sequencing Analysis Workflow 🧬

[![R](https://img.shields.io/badge/R-4.1.0-blue.svg)](https://www.r-project.org/)
[![Seurat](https://img.shields.io/badge/Seurat-4.0.0-green.svg)](https://satijalab.org/seurat/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

## 📋 Overview

This repository contains a comprehensive single-cell RNA sequencing (scRNA-seq) analysis  Which i have used to process the NYU samples for  analyzing pancreatic cancer samples. The analysis is based on data from the Gene Expression Omnibus (GEO) database (accession number GSE205013), which includes 27 samples (P01-P27). The analysis integrates multiple samples, performs quality control, removes batch effects, and identifies cell types using marker genes.

This repository contains a comprehensive single-cell RNA sequencing (scRNA-seq) analysis is designed for processing and analyzing pancreatic cancer samples. This analysis integrates multiple samples, performs quality control, removes batch effects, and identifies cell types using marker genes.

## 🔍 Features

- **Data Processing**
  - Automated processing of multiple samples (P01-P27)
  - Creation of Seurat objects from raw count matrices
  - Quality control and filtering
  - Removal of erythroid cells
  - Doublet detection and removal using scDblFinder

- **Analysis Pipeline**
  - Data normalization and feature selection
  - Dimensionality reduction (PCA and UMAP)
  - Cell clustering with multiple resolution testing
  - Integration of multiple samples
  - Comprehensive cell type annotation
  - Differential expression analysis

- **Visualization**
  - UMAP visualizations
  - Feature plots for marker genes
  - Violin plots for gene expression
  - Heatmaps of top markers
  - Cluster proportion analysis

## 🛠️ Requirements

### R Packages
```R
library(Seurat)
library(tidyverse)
library(SingleCellExperiment)
library(scDblFinder)
```

### Data Structure
The Analysis data general structure format:
```
data_dir/
├── GSM[number]_P[XX]_matrix.mtx.gz
├── GSM[number]_P[XX]_features.tsv.gz
└── GSM[number]_P[XX]_barcodes.tsv.gz
```

## 🚀 Usage

1. **Data Processing**
```R
# Create Seurat objects
samples <- sprintf("P%02d", 1:27)
seurat_list <- list()
for(i in 1:length(samples)) {
    seurat_list[[i]] <- create_seurat(matrix_file, features_file, barcodes_file, sample)
}
```

2. **Quality Control**
```R
# Apply QC filters
merged_seurat <- subset(merged_seurat, 
                       subset = nFeature_RNA > 500 & 
                               nCount_RNA > 1500 & 
                               percent.mt < 15)
```

3. **Integration and Clustering**
```R
# Find integration anchors
anchors <- FindIntegrationAnchors(object.list = seurat_list, 
                                 anchor.features = features,
                                 dims = 1:30)

# Create integrated data
integrated_seurat <- IntegrateData(anchorset = anchors, dims = 1:30)
```
# PDAC Data Summary

This dataset contains **32,736 features** across **141,186 single cells**, processed . Below is the breakdown of samples and their respective single-cell counts.

## Seurat Object Details
- **Features**: 32,736
- **Single Cells**: 141,186
- **Assays**: 2
  - **RNA**: 30,736 features
  - **Integrated**: 1 other assay
- **Dimensional Reductions**: PCA, UMAP

## Single Cell Distribution

| Samples  | Single Cells |
|:---------|-------------:|
| P01      |         4521 |
| P02      |         3291 |
| P03      |        14157 |
| P04      |         2163 |
| P05      |         8995 |
| P06      |         1102 |
| P07      |         3544 |
| P08      |         4195 |
| P09      |         1210 |
| P10      |         3318 |
| P11      |         4242 |
| P12      |         2884 |
| P13      |         1920 |
| P14      |        11260 |
| P15      |         6472 |
| P16      |         2694 |
| P17      |         3804 |
| P18      |        11103 |
| P19      |         9370 |
| P20      |         5638 |
| P21      |         1913 |
| P22      |         4599 |
| P23      |         7118 |
| P24      |         2358 |
| P25      |         7655 |
| P26      |         3950 |
| P27      |         7710 |
| **Total**|     **141186** |

## 📊 Cell Type Markers

I have used the following comprehensive marker sets for cell type identification:

```R
cell_markers <- list(
- T/NK cells: CD3E, CD3D, CD3G, NKG7, KLRD1, NCAM1, XCL1, GNLY
- Epithelial cells: KRT19, EPCAM, KRT8, KRT18, TFF1, TFF2, TFF3
- Endothelial cells: VWF, PECAM1, CDH5, PLVAP, ENG
- Myeloid cells: CD68, LYZ, CD14, CD163, C1QA, C1QB, C1QC, SPP1
- CAFs: DCN, COL1A1, COL3A1, ACTA2, MMP11, C3, C7, CFD, PTGDS
- B/Plasma cells: CD79A, CD79B, MS4A1, SDC1, MZB1, JCHAIN
- Mast cells: KIT, TPSAB1, CPA3, HPGDS
- Proliferating cells: MKI67, TOP2A, PCNA
# ... 
)
```

## 📈 Output Files

I have generated several output files:
- `PDAC_cell_type_annotation.pdf`: Comprehensive visualization report
- `cluster_markers_statistics.csv`: Detailed marker statistics
- `pancreas_treated_vs_naive_DE.csv`: Differential expression analysis results
- `annotated_seurat.rds`: Final annotated Seurat object
- `integrated_seurat.rds`: Integrated Seurat object

## 🤝 Contributing

Feel free to submit issues, fork the repository, and create pull requests for any improvements.

## 📝 License

This project is licensed under the MIT License - Connect with me (atheeq.khaliq@gmail.com)

## 📚 Citation

please cite:
```
Werba, G., Weissinger, D., Kawaler, E.A. et al. Single-cell RNA sequencing reveals the effects of chemotherapy on human pancreatic adenocarcinoma and its tumor microenvironment. Nat Commun 14, 797 (2023). https://doi.org/10.1038/s41467-023-36296-4
```

## 📧 Contact

For questions or collaborations, please contact:
- Email: akhaliq@iu.edu

Methods:

# Methods

## Single-Cell RNA Sequencing Data Processing and Analysis

### Data Source
Single-cell RNA sequencing data was obtained from the Gene Expression Omnibus (GEO) database under accession number GSE205013. Raw data files for 27 samples (P01-P27) were downloaded and processed according to the following protocol.

### Data Processing and Quality Control
Single-cell RNA sequencing data processing was performed using R (version 4.1.0) with the Seurat package (version 4.0.0). Raw count matrices were processed using the CreateSeuratObject function with initial filtering criteria of minimum 3 cells per gene and minimum 200 features per cell. Quality control metrics were calculated for each cell, including the percentage of mitochondrial genes (pattern "^MT-"). Cells were filtered based on the following criteria: number of features > 500, total counts > 1,500, and mitochondrial percentage < 15%. 

Potential contaminating erythroid cells were removed by calculating the percentage of expression for erythroid-specific genes (HBA1, HBA2, HBB, HBM, ALAS2) and excluding cells with erythroid gene expression > 1% of total counts. Doublet detection was performed using the scDblFinder package (version 1.8.0), taking into account sample-specific information.

### Data Integration and Dimensionality Reduction
Data normalization was performed using the LogNormalize method in Seurat. Variable features were identified using the FindVariableFeatures function with the "vst" method, selecting the top 2,000 variable genes. Data integration across samples was performed using Seurat's integration workflow. Integration anchors were identified using FindIntegrationAnchors (dims = 1:30), and data integration was performed using IntegrateData. 

Principal Component Analysis (PCA) was performed on the integrated data using RunPCA. The first 30 principal components were used for downstream analysis. Uniform Manifold Approximation and Projection (UMAP) was performed using RunUMAP with parameters: dims = 1:30, n.neighbors = 30, and min.dist = 0.3.

### Clustering and Cell Type Annotation
Graph-based clustering was performed using FindNeighbors (dims = 1:30) followed by FindClusters. Multiple resolution parameters (0.1-10.0) were tested, with a final resolution of 0.8 selected based on biological relevance. The Leiden algorithm (algorithm = 4) was used for community detection.

Cell type annotation was performed using established marker genes for major cell populations:
- T/NK cells: CD3E, CD3D, CD3G, NKG7, KLRD1, NCAM1, XCL1, GNLY
- Epithelial cells: KRT19, EPCAM, KRT8, KRT18, TFF1, TFF2, TFF3
- Endothelial cells: VWF, PECAM1, CDH5, PLVAP, ENG
- Myeloid cells: CD68, LYZ, CD14, CD163, C1QA, C1QB, C1QC, SPP1
- CAFs: DCN, COL1A1, COL3A1, ACTA2, MMP11, C3, C7, CFD, PTGDS
- B/Plasma cells: CD79A, CD79B, MS4A1, SDC1, MZB1, JCHAIN
- Mast cells: KIT, TPSAB1, CPA3, HPGDS
- Proliferating cells: MKI67, TOP2A, PCNA

Module scores for each cell type were calculated using AddModuleScore. Additional subtype-specific markers were used for detailed annotation of T/NK cells, myeloid cells, and CAFs.

### Differential Expression Analysis
Differential expression analysis between treatment groups (Naïve vs Treated) was performed using Seurat's FindMarkers function with the Wilcoxon rank-sum test. Analysis was restricted to pancreatic samples. Genes were considered significantly differentially expressed if they met the following criteria: adjusted p-value < 0.05 and |log2 fold change| > 0.25. A minimum detection threshold of 10% (min.pct = 0.1) was applied to ensure robust differential expression analysis.

### Software and Packages
I have utilized the following R packages:
- Seurat (v4.4.0)
- tidyverse (v2.0.0)
- SingleCellExperiment (v1.28.1)
- scDblFinder (v1.20.0)
- ggplot2 (v3.5.1)

### Data Visualization
Visualization of results included UMAP plots for dimensionality reduction, feature plots for marker gene expression, violin plots for gene expression distribution, dot plots for marker gene expression across clusters, and heatmaps for differentially expressed genes. All visualizations were generated using ggplot2 and Seurat's built-in plotting functions.

### Data Availability
The final processed data has been saved in multiple formats:
- Annotated Seurat object (annotated_seurat.rds)
- Integrated Seurat object (integrated_seurat.rds)
- Differential expression results (pancreas_treated_vs_naive_DE.csv)
- Cluster marker statistics (cluster_markers_statistics.csv)
- Comprehensive visualization report (PDAC_cell_type_annotation.pdf)

### Reference
Werba, G., Weissinger, D., Kawaler, E.A. et al. Single-cell RNA sequencing reveals the effects of chemotherapy on human pancreatic adenocarcinoma and its tumor microenvironment. Nat Commun 14, 797 (2023). https://doi.org/10.1038/s41467-023-36296-4

