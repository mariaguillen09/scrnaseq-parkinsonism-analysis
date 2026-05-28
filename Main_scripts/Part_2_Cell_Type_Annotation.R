# ==============================================================================
# Script: Part_2_Cell_Type_Annotation.R
# Project: Single-cell transcriptomic profiling of peripheral immune cells in 
#          genetically stratified Parkinson’s disease and progressive supranuclear 
#          palsy
# Description: This script performs automatic and manual cell type annotation,
#              quality filtering based on prediction confidence, and validation
#              using canonical marker genes.
# ==============================================================================

# Load Required Libraries
# scRNA-seq Analysis
library(Seurat)

# Automated annotation tools
library(SingleR)      # Automatic cell type annotation
library(celldex)      # Reference datasets for SingleR
library(Azimuth)      # Azimuth automatic annotation (PBMC reference)
library(Signac)       # Companion package for Azimuth
library(SeuratData)   # Data resources for Seurat

# Visualization
library(pheatmap)
library(ggplot2)
library(patchwork)

# Data handling and Utilities
library(openxlsx)
library(tidyverse)

### 1. Load Integrated Object from Part 1.
obj_final <- readRDS("Results/Integrated_Seurat_Object.rds")

# Join layers for annotation
obj_final <- JoinLayers(obj_final)

### 2. Automatic Cell Type Annotation

## 2.1 SingleR Annotation
# SingleR uses a reference dataset (Monaco Immune Data) to assign cell
# types to each cell based on similarity of gene expression patterns.

# Load reference dataset
monaco <- celldex::MonacoImmuneData()

# Run SingleR annotation
pred.monaco <- SingleR(
  test = obj_final@assays$RNA$data,  
  ref = monaco,                       
  labels = monaco$label.main          
)

# Add predictions to Seurat metadata
obj_final$singleR_annotation <- pred.monaco$labels

cat(sprintf("  Cell types identified: %d\n", nlevels(factor(obj_final$singleR_annotation))))

# Show distribution
cat("\nCell type distribution (SingleR):\n")
print(table(obj_final$singleR_annotation))

# Visualize SingleR annotations
singleR_plot <- DimPlot(
  obj_final, 
  group.by = "singleR_annotation", 
  label = TRUE,
  label.size = 3
) + 
  ggtitle("SingleR Automatic Annotation")

ggsave(
  filename = "Figures/09_SingleR_Annotation.png",
  plot = singleR_plot,
  width = 12,
  height = 8,
  dpi = 300
)

## 2.2 Azimuth Annotation
# Azimuth uses a pre-trained model on human PBMC data to automatically
# assign cell types and also provides confidence scores for each
# prediction.

# Run Azimuth annotation
obj_final <- RunAzimuth(obj_final, reference = "pbmcref")

# Azimuth provides multiple levels of cell type resolution
# l1 = broad cell types
cat("\nAzimuth cell type predictions (level 1):\n")
print(table(obj_final$predicted.celltype.l1))

# Check prediction confidence
cat("\nAzimuth prediction confidence (l1):\n")
cat(sprintf("  Mean confidence: %.3f\n", mean(obj_final$predicted.celltype.l1.score)))
cat(sprintf("  Cells with confidence > 0.5: %d (%.1f%%)\n",
            sum(obj_final$predicted.celltype.l1.score > 0.5),
            (sum(obj_final$predicted.celltype.l1.score > 0.5) / ncol(obj_final)) * 100))

# Visualize Azimuth annotations
azimuth_plot <- DimPlot(
  obj_final, 
  group.by = "predicted.celltype.l1", 
  label = TRUE,
  repel = TRUE
) + 
  ggtitle("Azimuth Automatic Annotation")

ggsave(
  filename = "Figures/10_Azimuth_Annotation.png",
  plot = azimuth_plot,
  width = 12,
  height = 8,
  dpi = 300
)

## 2.3 Quality Filter Based on Azimuth Confidence
# Cells with low prediction confidence are likely doublets or cell types
# not well-represented in the reference.

# Keep only cells with confidence prediction score > 0.5
cells_to_keep <- obj_final$predicted.celltype.l1.score > 0.5
n_removed <- sum(!cells_to_keep)

obj_final <- subset(obj_final, cells = names(which(cells_to_keep)))

cat(sprintf("  Cells removed: %d (likely doublets/uncertain types)\n", n_removed))
cat(sprintf("  Cells remaining: %d\n", ncol(obj_final)))

### 3. Identify Cluster-Specific Marker Genes
# To manually annotate clusters, we identify the genes that are most
# highly expressed in each cluster. These marker genes help us understand
#the identity of each cluster.

## 3.1 Find All Cluster Markers

# Find all markers for each cluster
# Only positive markers (genes up-regulated in the cluster)
all_markers <- FindAllMarkers(
  obj_final,
  only.pos = TRUE,           # Only upregulated genes
  min.pct = 0.25,            # Gene must be in 25% of cluster cells
  logfc.threshold = 0.25,    # Minimum log2 fold change
  verbose = FALSE
)

# Get top 10 markers per cluster for manual review
top10_markers <- all_markers %>%
  group_by(cluster) %>%
  slice_max(n = 10, order_by = avg_log2FC)

write.xlsx(top10_markers, "Results/02_Top10_ClusterMarkers.xlsx")

## 3.2 Manual Cell Type Assignment
# Based on the marker genes and automatic annotation, we manually assign
# cell types to clusters.

# Manual mapping of clusters to cell types
# This should be adjusted based on your marker gene analysis above
cluster_to_celltype <- c(
  "0" = "T_CD4_cells",
  "1" = "T_CD4_cells",
  "2" = "Monocytes",
  "3" = "T_CD8_cells",
  "4" = "NK_cells",
  "5" = "T_CD8_cells",
  "6" = "Monocytes",
  "7" = "T_CD8_cells",
  "8" = "B_cells",
  "9" = "B_cells",
  "10" = "Platelets",
  "11" = "Platelets",
  "12" = "Dendritic_cells",
  "13" = "Dendritic_cells",
  "14" = "proliferative_cells",
  "15" = "Basophils",
  "16" = "HSC",
  "17" = "Neutrophils",
  "18" = "B_cells"
)

# Map cluster identities to cell types
obj_final$final_annotation <- plyr::mapvalues(
  x = obj_final$seurat_clusters,
  from = names(cluster_to_celltype),
  to = cluster_to_celltype
)

obj_final$final_annotation <- factor(
  obj_final$final_annotation,
  levels = c(
    "T_CD4_cells", "T_CD8_cells", "NK_cells", "B_cells",
    "Monocytes", "Dendritic_cells", "HSC", "proliferative_cells",
    "Platelets", "Neutrophils", "Basophils"
  )
)

cat("\nCell type distribution:\n")
print(table(obj_final$final_annotation))

### 4. Validation with Marker Gene Expression
# We create a dot plot showing the expression of canonical marker genes
# across cell types to validate cluster annotations.

## 4.1 Define Canonical Markers

# Define canonical markers for validation
# These are well-known genes that identify each cell type
marker_genes <- c(
  # T cell markers
  "CD4", "CCR4", "GATA3", "CD40LG",        # CD4+ T cells
  "CD8A", "GZMK", "CD8B",                   # CD8+ T cells
  # NK cells
  "NCAM1", "SPON2", "CLIC3",
  # B cells
  "CD19", "IGHD", "CD22",
  # Monocytes
  "CD14", "VCAN", "CDKN1C",
  # Dendritic cells
  "CLEC4C", "FCER1A", "CD1C",
  # Other
  "CD34", "MKI67", "PPBP", "FCGR3B", "HDC"
)

## 4.2 Create Validation Dot Plot
# Create dot plot showing marker expression per cell type

validation_plot <- DotPlot(
  obj_final,
  features = marker_genes,
  group.by = "final_annotation",
  cols = c("lightgrey", "#D73027")
) +
  coord_flip() +
  theme_bw(base_size = 20) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 32, margin = margin(b = 25)),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, face = "bold", size = 25),
    axis.text.y = element_text(face = "italic", size = 25),
    axis.title.x = element_text(size = 25, face = "bold", margin = margin(t = 15)),
    axis.title.y = element_text(size = 25, face = "bold", margin = margin(r = 15)),
    legend.title = element_text(size = 25, face = "bold"),
    legend.text = element_text(size = 25),
    legend.position = "right",
    panel.grid.major = element_line(color = "grey95", linewidth = 0.2)
  ) +
  labs(
    title = "Validation of Canonical Cell Type Markers",
    x = "Marker Genes",
    y = "Cell Type",
    color = "Mean Expression",
    size = "% Cells Expressing"
  )


ggsave(
  filename = "Figures/11_MarkerValidation_DotPlot.png",
  plot = validation_plot,
  width = 13,
  height = 18,
  dpi = 300
)

## 4.3 UMAP Visualization with Final Annotations
umap_clusters <- DimPlot(
  obj_final,
  reduction = "umap",
  group.by = "seurat_clusters",
  label = TRUE, label.size = 7
) +
  ggtitle("UMAP - Seurat Clusters") +
  theme(legend.position = "right") + 
  theme(axis.text = element_text(size = 20),
        axis.title = element_text(size = 20),
        plot.title = element_text(size = 25, face = "bold"),
        legend.text = element_text(size = 20),
        legend.key.height = unit(1, "cm")) +
  guides(color = guide_legend(override.aes = list(size = 5)))

umap_celltypes <- DimPlot(
  obj_final,
  reduction = "umap",
  group.by = "final_annotation",
  label = TRUE,label.size = 7,repel = TRUE
) +
  ggtitle("UMAP - Cell Type Annotation") +
  theme(legend.position = "right") + 
  theme(axis.text = element_text(size = 20),
        axis.title = element_text(size = 20),
        plot.title = element_text(size = 25, face = "bold"),
        legend.text = element_text(size = 20),
        legend.key.height = unit(1, "cm")) +
  guides(color = guide_legend(override.aes = list(size = 5)))

# Combine into single figure
combined_umap <- (umap_clusters / umap_celltypes)

ggsave(
  filename = "Figures/12_Combined_UMAP_Annotations.png",
  plot = combined_umap,
  width = 12,
  height = 18,
  dpi = 300
)

### 5. Create Diagnostic Comparison Groups
# For downstream analysis, we create a unified diagnostic grouping
# variable.

# Create comparison groups
obj_final$comparison_group <- as.character(obj_final$diagnosis)

# Standardize naming
obj_final$comparison_group[obj_final$comparison_group %in% 
                             c("PD_LRRK2", "PD_GBA", "PD-LRRK2", 
                               "PD-GBA")] <- "Genetic-PD"
obj_final$comparison_group[obj_final$comparison_group %in% 
                             c("sporadic_PD", 
                               "sporadic-PD")] <- "sporadic-PD"

# Convert to factor with logical order
obj_final$comparison_group <- factor(
  obj_final$comparison_group,
  levels = c("HC", "sporadic-PD", "Genetic-PD", "PSP")
)

### 7. Save Annotated Object
# Save the final annotated object for downstream analysis.
saveRDS(obj_final, file = "Results/Final_Annotated_Seurat_Object.rds")


