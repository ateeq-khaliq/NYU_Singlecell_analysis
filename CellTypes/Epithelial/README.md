# GSEA Analysis Pipeline

## Overview
This repository details about how i performed Gene Set Enrichment Analysis (GSEA) on differential expression data comparing treated vs. naive conditions. The pipeline includes differential expression analysis using Seurat, followed by GSEA using the Molecular Signatures Database (MSigDB) Hallmark gene sets, and generates publication-ready visualizations.

## Dependencies
The following R packages are required:
- Seurat (for differential expression analysis)
- clusterProfiler (for GSEA)
- msigdbr (for accessing MSigDB gene sets)
- dplyr (for data manipulation)
- ggplot2 (for visualization)

You can install these packages using:
```R
install.packages(c("Seurat", "clusterProfiler", "msigdbr", "dplyr", "ggplot2"))
```

## Workflow
1. **Differential Expression Analysis**
   - Uses Seurat's `FindMarkers` function
   - Compares "Treated" vs "Naive" conditions
   - Implements Wilcoxon rank sum test

2. **GSEA Preparation**
   - Calculates ranking metric using adjusted p-values and log2 fold changes
   - Sorts genes by ranking metric
   - Retrieves Hallmark gene sets for Homo sapiens

3. **GSEA Implementation**
   - Performs GSEA with following parameters:
     - Minimum gene set size: 10
     - Maximum gene set size: 500
     - P-value cutoff: 1 (all results shown)

4. **Visualization**
   - Creates a horizontal bar plot showing:
     - Normalized Enrichment Scores (NES)
     - Color-coded by enrichment direction (red for Treated, blue for Naive)
     - Statistical significance indicators (* p<0.05, ** p<0.01, *** p<0.001)

## Output
- Generates a PDF file named "previous_hallmark_pathways.pdf"
- Plot dimensions: 14x16 inches
- Includes all Hallmark pathways with their enrichment scores


## Visualization Features
- Horizontal bar plot
- Pathways ordered by NES
- Clear significance indicators
- Clean, minimal theme
- Dashed line at zero for reference
- Legend at bottom
- Customized text sizes for readability

## METHODS


Gene Set Enrichment Analysis of Differential Expression Data

Differential Expression Analysis
Differential expression analysis between treated and naive conditions was performed using Seurat's FindMarkers function implementing the Wilcoxon rank sum test. The analysis was conducted on a cell-type specific basis, focusing on the epithelial cell population. Genes with adjusted p-values < 0.05 were considered significantly differentially expressed.

Gene Set Enrichment Analysis
Gene Set Enrichment Analysis (GSEA) was performed using the clusterProfiler package in R. A ranking metric was calculated for each gene using the formula:
ranking_metric = -log10(adjusted p-value) × sign(log2 fold change)

This metric incorporated both statistical significance and the direction of expression changes. Genes were ranked in descending order based on this metric. The Molecular Signatures Database (MSigDB) Hallmark gene set collection (Liberzon et al., 2015) for Homo sapiens was obtained using the msigdbr package. GSEA was performed with the following parameters: minimum gene set size = 10, maximum gene set size = 500, and p-value cutoff = 1.0 to obtain a complete ranking of pathways.

Statistical Analysis
Multiple testing correction was performed using the Benjamini-Hochberg procedure to control the false discovery rate. Pathways with adjusted p-values < 0.05 were considered significantly enriched. Normalized Enrichment Scores (NES) were calculated to account for differences in gene set sizes and correlations between gene sets.

Visualization
Results were visualized using ggplot2 in R. Enrichment results were displayed as horizontal bar plots with the following features:
- Normalized Enrichment Scores (NES) represented by bar length
- Directionality of enrichment indicated by color (red for treated condition, blue for naive condition)
- Statistical significance denoted by asterisks (* p<0.05, ** p<0.01, *** p<0.001)
- Pathways ordered by NES magnitude
- Dashed reference line at NES = 0



References
1. Stuart T, Butler A, Hoffman P, et al. Comprehensive Integration of Single-Cell Data. Cell. 2019
2. Wu T, Hu E, Xu S, et al. clusterProfiler 4.0: A universal enrichment tool for interpreting omics data. Innovation. 2021
3. Liberzon A, Birger C, Thorvaldsdóttir H, et al. The Molecular Signatures Database Hallmark Gene Set Collection. Cell Systems. 2015
4. Wickham H. ggplot2: Elegant Graphics for Data Analysis. Springer-Verlag New York, 2016
