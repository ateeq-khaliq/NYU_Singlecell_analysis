# Load required libraries
library(Seurat)
library(dplyr)
library(ggplot2)

# Load the CAF subset Seurat object
caf_obj <- readRDS("/Users/akhaliq/Desktop/sc_conformation_data/nyu_data/celltypes/pdac_pancreas_cafs.rds")

# Load required libraries
library(Seurat)
library(dplyr)
library(ggplot2)

# Load the CAF subset Seurat object
#caf_obj <- readRDS("pdac_pancreas_cafs.rds")

# Make sure we're using RNA assay for marker analysis
DefaultAssay(caf_obj) <- "RNA"

# Define core marker genes for each subtype
caf_markers_original <- list(
  myCAF = c("MYL9", "ACTA2", "TAGLN", "TPM1", "TPM2", "POSTN", "CTHRC1", "BGN", "COL10A1","TPM2"),
  iCAF = c("CXCL12", "PDGFRA", "CFD", "CXCL14", "FBLN1", "DPT", "APOE", "PTGDS", "CCL2", "IL6"),
  apCAF = c( "SLPI", "CD74", "HLA-DRA", "HLA-DPA1", "HLA-DPB1", "HLA-DQB1")
)

# Get available features in the dataset
available_features <- rownames(caf_obj)

# Function to filter markers and print availability status
filter_markers <- function(marker_list, available_features) {
  present_markers <- intersect(marker_list, available_features)
  missing_markers <- setdiff(marker_list, available_features)
  
  if(length(missing_markers) > 0) {
    print(paste("Missing markers for", deparse(substitute(marker_list)), ":"))
    print(missing_markers)
  }
  
  return(present_markers)
}

# Filter markers to only include available genes
caf_markers <- list(
  myCAF = filter_markers(caf_markers_original$myCAF, available_features),
  iCAF = filter_markers(caf_markers_original$iCAF, available_features),
  apCAF = filter_markers(caf_markers_original$apCAF, available_features)
)

# Print number of markers used for each type
print("\nNumber of markers used for each CAF subtype:")
print(sapply(caf_markers, length))

# Calculate module scores for each CAF subtype
print("\nCalculating module scores...")
caf_obj <- AddModuleScore(caf_obj, 
                         features = caf_markers$myCAF,
                         name = "myCAF_score",
                         ctrl = 100)

caf_obj <- AddModuleScore(caf_obj, 
                         features = caf_markers$iCAF,
                         name = "iCAF_score",
                         ctrl = 100)

caf_obj <- AddModuleScore(caf_obj, 
                         features = caf_markers$apCAF,
                         name = "apCAF_score",
                         ctrl = 100)

# Simple classification based on highest score
scores <- cbind(caf_obj$myCAF_score1,
                caf_obj$iCAF_score2,
                caf_obj$apCAF_score3)
colnames(scores) <- c("myCAF", "iCAF", "apCAF")

caf_obj$caf_subtype <- colnames(scores)[max.col(scores)]

# Create visualizations
print("\nCreating visualizations...")

# UMAP with CAF subtypes
p1 <- DimPlot(caf_obj, 
              reduction = "umap",
              group.by = "caf_subtype",
              label = TRUE,
              label.size = 4) +
  ggtitle("CAF Subtypes") +
  theme_minimal()

ggsave("caf_subtypes_umap.pdf", p1, width = 10, height = 8)

# Feature plots for module scores
p2 <- FeaturePlot(caf_obj, 
                  features = c("myCAF_score1", "iCAF_score2", "apCAF_score3"),
                  ncol = 3,
                  order = TRUE)

ggsave("caf_subtype_scores.pdf", p2, width = 15, height = 5)

# Feature plots for individual key markers
key_markers <- list(
  myCAF = c("ACTA2", "TAGLN", "TPM1"),
  iCAF = c("CXCL12", "IL6", "PTGDS"),
  apCAF = c("KRT19", "CD74", "HLA-DRA")
)

# Plot key markers
for(caf_type in names(key_markers)) {
  markers <- intersect(key_markers[[caf_type]], available_features)
  if(length(markers) > 0) {
    p <- FeaturePlot(caf_obj, 
                     features = markers,
                     ncol = 3,
                     order = TRUE)
    ggsave(paste0("caf_", tolower(caf_type), "_markers.pdf"), p, width = 15, height = 5)
  }
}

# Create violin plots for key markers
all_key_markers <- unlist(key_markers)
all_key_markers <- intersect(all_key_markers, available_features)

p3 <- VlnPlot(caf_obj, 
              features = all_key_markers,
              group.by = "caf_subtype",
              ncol = 3,
              pt.size = 0)

ggsave("caf_key_markers_violin.pdf", p3, width = 15, height = 15)

# Create heatmap
avg_expr <- AverageExpression(caf_obj, 
                             features = unique(unlist(caf_markers)),
                             group.by = "caf_subtype")

# Scale the expression data
scaled_expr <- t(scale(t(avg_expr$RNA)))

# Create heatmap
library(pheatmap)
pdf("caf_subtype_heatmap.pdf", width = 12, height = 20)
pheatmap(scaled_expr,
         cluster_rows = TRUE,
         cluster_cols = FALSE,
         fontsize_row = 8,
         main = "CAF Subtype Marker Expression")
dev.off()

# Save results
saveRDS(caf_obj, "pdac_pancreas_cafs_classified.rds")

# Export cell assignments to CSV
cell_assignments <- data.frame(
  Cell_ID = colnames(caf_obj),
  Original_Cluster = caf_obj$seurat_clusters,
  CAF_Subtype = caf_obj$caf_subtype,
  myCAF_Score = caf_obj$myCAF_score1,
  iCAF_Score = caf_obj$iCAF_score2,
  apCAF_Score = caf_obj$apCAF_score3
)

write.csv(cell_assignments, "caf_subtype_assignments.csv", row.names = FALSE)

# Print summary statistics
print("\nCAF subtype distribution:")
print(table(caf_obj$caf_subtype))
print("\nPercentage of each subtype:")
print(round(prop.table(table(caf_obj$caf_subtype)) * 100, 2))

# Save marker gene lists actually used
write.csv(data.frame(
  Gene = unlist(caf_markers),
  Type = rep(names(caf_markers), sapply(caf_markers, length))
), "used_marker_genes.csv", row.names = FALSE)

