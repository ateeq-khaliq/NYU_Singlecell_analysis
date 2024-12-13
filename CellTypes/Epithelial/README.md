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
