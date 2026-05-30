
# ==============================================================================
# Script: Part_1_Preprocessing_QC_Integration.R
# Project: Single-cell transcriptomic profiling of peripheral immune cells in 
#          genetically stratified Parkinson’s disease and progressive supranuclear 
#          palsy
# Description: This script performs raw data loading, quality control (QC), 
#              doublet removal, normalization, dimensionality reduction, 
#              and benchmarks batch effect correction methods (Harmony, RPCA, 
#              FastMNN)
#              using LISI metrics to select the optimal integration strategy.
# ==============================================================================

# Create output directories
dir.create("Figures")
dir.create("Results")

# Load Required Libraries
# packages for single-cell analysis
library(Seurat)          
library(tidyverse)       
library(ggplot2)         
library(harmony)        
library(SeuratWrappers) 
library(batchelor)     
library(scDblFinder)  
library(lisi)      
library(future)  
library(openxlsx)        
options(future.globals.maxSize = 100 * 1024^3)

### 1. Data Loading and Initial Setup
# In this section, we load the raw count matrices from 10X Genomics format and 
# create individual Seurat objects for each sample.

## 1.1 Define Sample Metadata
# First, we organize information about all samples we're analyzing. This includes 
# the diagnosis, file paths, pool assignment (batch), and unique identifiers.

# Base path to sequencing data
base_path <- "./data_tfm/"

# Create metadata table with sample information
# This table links each sample to its diagnosis, file location, and experimental 
# batch
samples_meta <- data.frame(
  diagnosis = c(
    # Pool A
    "sporadic_PD", "PD_LRRK2", "PSP", "HC",           
    # Pool B
    "sporadic_PD", "PSP", "PD_GBA", "HC",
    # Pool C
    "PD_LRRK2", "PSP", "sporadic_PD", "HC",
    # Pool D
    "PD_GBA", "PSP", "HC", "sporadic_PD"
  ),
  
  file_path = paste0(base_path, c(
    # Pool A
    "poolA/16646", "poolA/16093", "poolA/16202", "poolA/17646",
    # Pool B
    "poolB/16131", "poolB/16740", "poolB/17056", "poolB/17878",
    # Pool C
    "poolC/16686", "poolC/17907", "poolC/18427", "poolC/18437",
    # Pool D
    "poolD/14127", "poolD/16414", "poolD/16727", "poolD/18050"
  )),
  
  pool = c(
    rep("Pool_A", 4), rep("Pool_B", 4), 
    rep("Pool_C", 4), rep("Pool_D", 4)
  ),
  
  id = c(
    "sPD_A", "LRRK2_A", "PSP_A", "HC_A", 
    "sPD_B", "PSP_B", "GBA_B", "HC_B",
    "LRRK2_C", "PSP_C", "sPD_C", "HC_C",
    "GBA_D", "PSP_D", "HC_D", "sPD_D"
  )
)

## 1.2 Load Count Matrices
# Now we read the raw count matrices and create Seurat objects for each sample.

# Initialize empty list to store Seurat objects
seurat_pre_qc <- list()

# Loop through each sample and load its data
for (i in seq_along(samples_meta$id)) {
  id <- samples_meta$id[i]
  
  # Read 10X format count matrix
  data <- Read10X(samples_meta$file_path[i])
  
  # Create Seurat object
  # min.cells = 3 means a gene must be detected in at least 3 cells to be kept
  obj <- CreateSeuratObject(
    counts = data, 
    project = samples_meta$diagnosis[i],
    min.cells = 3
  )
  
  # Add metadata to track experimental information
  obj$pool <- samples_meta$pool[i]           # Which batch/pool sample came from
  obj$id_tecnico <- id                       # Sample identifier
  obj$diagnosis <- samples_meta$diagnosis[i] # Disease status
  
  # Store in list
  seurat_pre_qc[[id]] <- obj
}

# Summary statistics
total_cells_initial <- sum(sapply(seurat_pre_qc, ncol))
cat(sprintf("  Total cells pre-QC: %d\n", total_cells_initial))

### 2. Quality Control and Doublet Detection

# Quality control is crucial for removing low-quality cells that could bias 
# downstream analysis. We'll calculate QC metrics, visualize them, apply 
# filtering thresholds, and remove doublets (cells where two cells were 
# captured in one droplet).

## 2.1 Calculate QC Metrics
# For each cell, we calculate nFeature_RNA, nCount_RNA, percent.mt and percent.hb.

# Loop through samples and add QC metrics
for (id in names(seurat_pre_qc)) {
  
  # Calculate mitochondrial percentage (MT genes are markers of cell stress)
  seurat_pre_qc[[id]][["percent.mt"]] <- PercentageFeatureSet(
    seurat_pre_qc[[id]], 
    pattern = "^MT-"
  )
  
  # Calculate hemoglobin percentage (RBC contamination marker)
  # Only include hemoglobin genes that are actually present in the data
  hb_genes <- c("HBA1", "HBA2", "HBB", "HBD", "HBE1", "HBG1", "HBG2", "HBM", 
                "HBQ1", "HBZ")
  hb_genes <- intersect(hb_genes, rownames(seurat_pre_qc[[id]]))
  seurat_pre_qc[[id]][["percent.hb"]] <- PercentageFeatureSet(
    seurat_pre_qc[[id]], 
    features = hb_genes
  )
}

## 2.2 Merge and Visualize Pre-QC Distribution
# To visualize QC metrics across all samples together, we merge them and decide 
# filtering thresholds.

# Merge individual Seurat Objects
obj_merged <- merge(
  seurat_pre_qc[[1]], 
  y = seurat_pre_qc[-1], 
  add.cell.ids = names(seurat_pre_qc)
)

# Create violin plot showing distributions grouped by diagnosis
pre_qc_plot <- VlnPlot(
  obj_merged, 
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.hb"), 
  ncol = 4, 
  group.by = "diagnosis", 
  pt.size = 0
) & 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Save the pre-QC visualization
ggsave(
  filename = "Figures/01_ViolinPlots_PreQC.png", 
  plot = pre_qc_plot, 
  width = 20,         
  height = 8,          
  dpi = 300
)

cat(sprintf("  Total cells before filtering: %d\n", ncol(obj_merged)))
cat("\nCells per diagnosis:\n")
print(table(obj_merged$diagnosis))

## 2.3 Apply QC Filters
# Based on the distributions visualized above, we apply filtering thresholds to 
# remove low-quality cells and doublets.

# Apply QC filtering with carefully chosen thresholds
obj_filtered <- subset(
  obj_merged, 
  subset = nFeature_RNA > 200 &       # Minimum genes per cell
    nFeature_RNA < 5500 &       # Maximum genes per cell 
    nCount_RNA < 15000 &        # Maximum UMI count
    percent.mt < 5 &            # Maximum mitochondrial content
    percent.hb < 5              # Maximum hemoglobin content
)

# Visualize post-QC distributions
post_qc_plot <- VlnPlot(
  obj_filtered, 
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.hb"), 
  ncol = 4, 
  group.by = "diagnosis", 
  pt.size = 0
) & 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(
  filename = "Figures/02_ViolinPlots_PostQC.png", 
  plot = post_qc_plot, 
  width = 20,         
  height = 8,          
  dpi = 300   
)

cat(sprintf("  Total cells after QC filtering: %d\n", ncol(obj_filtered)))
cat(sprintf("  Cells removed: %d (%.1f%% of original)\n", 
            ncol(obj_merged) - ncol(obj_filtered),
            (1 - ncol(obj_filtered)/ncol(obj_merged)) * 100))
cat("\nCells per diagnosis after filtering:\n")
print(table(obj_filtered$diagnosis))

## 2.4 Detect and Remove Doublets
# Doublets are two cells captured in the same droplet. They can bias downstream 
# analysis. We use scDblFinder to detect and remove them.

# Loop through original samples (not the merged object) for accurate doublet detection
for (id in names(seurat_pre_qc)) {
  
  # Convert Seurat object to SingleCellExperiment format (required by scDblFinder)
  sce <- as.SingleCellExperiment(seurat_pre_qc[[id]])
  
  # Run doublet detection
  sce <- scDblFinder(sce)
  
  seurat_pre_qc[[id]] <- RenameCells(seurat_pre_qc[[id]], 
  new.names = paste0(id, "_", colnames(seurat_pre_qc[[id]])))
  
  n_doublets <- sum(sce$scDblFinder.class == "doublet")
  cat(sprintf("  %s: %d doublets detected\n", 
              id, n_doublets))
  
  # Add doublet information back to Seurat object
  seurat_pre_qc[[id]]$doublet_class <- sce$scDblFinder.class
}

todas_las_etiquetas <- unlist(lapply(seurat_pre_qc, function(x) x$doublet_class), 
                              use.names = FALSE)
names(todas_las_etiquetas) <- unlist(lapply(seurat_pre_qc, function(x) colnames(x)), 
                                     use.names = FALSE)

obj_filtered$doublet_class <- todas_las_etiquetas[colnames(obj_filtered)]

# Filter out doublets from the QC-filtered object
obj_filtered <- subset(obj_filtered, subset = doublet_class == "singlet")

cat(sprintf("  Total cells after doublet removal: %d\n", ncol(obj_filtered)))
table(obj_filtered$diagnosis)
cat(sprintf("  Cumulative cells retained: %d (%.1f%% of original)\n", 
            ncol(obj_filtered),
            (ncol(obj_filtered)/total_cells_initial) * 100))

### 3. Normalization, Feature Selection and Dimensionality Reduction

# Now that we have high-quality cells, we prepare the data for analysis by: 
# 1. Normalizing counts to account for sequencing depth differences 
# 2. Finding highly variable genes (features with biological signal) 
# 3. Scaling the data 
# 4. Running PCA for dimensionality reduction

## 3.1 Prepare RNA Layers

# In Seurat v5, RNA data can be split into separate layers. We join them 
# and then split by pool for normalization, HVG selection and batch effect 
# correction.

# Join all RNA layers
obj_filtered[["RNA"]] <- JoinLayers(obj_filtered[["RNA"]])

# Split by pool to correct batch effects
obj_filtered[["RNA"]] <- split(
  obj_filtered[["RNA"]], 
  f = obj_filtered$pool
)

## 3.2 Normalize Data
# Log-normalization adjusts for differences in sequencing depth between cells. 
# Each cell's counts are normalized to the same total (scale factor = 10,000), 
# then log-transformed to stabilize variance.

# Apply log-normalization
# This makes counts more comparable across cells and samples
obj_filtered <- NormalizeData(
  obj_filtered, 
  normalization.method = "LogNormalize", 
  scale.factor = 10000, 
  verbose = FALSE
)

## 3.3 Find Highly Variable Genes
# Identification of the 2,000 genes with the highest biological variance across 
# cells. These genes contain most of the meaningful biological signal.

# Find highly variable genes
obj_filtered <- FindVariableFeatures(
  obj_filtered,
  selection.method = "vst",  # Variance-stabilizing transformation
  nfeatures = 2000,
  verbose = FALSE
)

# Create HVG visualization
hvg_plot <- VariableFeaturePlot(obj_filtered)
hvg_plot_labeled <- LabelPoints(
  plot = hvg_plot, 
  points = head(VariableFeatures(obj_filtered), 10),  # top 10
  repel = TRUE
)

# Save HVG plot
ggsave(
  filename = "Figures/03_HVG_Plot.png", 
  plot = hvg_plot_labeled, 
  width = 10, 
  height = 8, 
  dpi = 300
)

## 3.4 Scale Data and Run PCA
# Scaling centers each gene's expression around zero and scales to unit variance. 
# PCA then projects this high-dimensional data into a lower-dimensional space 
# while preserving the main sources of variance.

# Scale data
obj_filtered <- ScaleData(obj_filtered)

# PCA
obj_filtered <- RunPCA(obj_filtered)

# Optimal number of PCs based on variance explained
stdevs <- obj_filtered[["pca"]]@stdev
percent.var <- (stdevs^2 / sum(stdevs^2)) * 100
cum_var <- cumsum(percent.var)

# PCs needed for 90% of variance
pcs_90 <- which(cum_var >= 90)[1]

cat(sprintf("  Number of PCs for 90%% variance: %d\n", pcs_90))

# Elbow plot to visualize variance explained
elbow_plot <- ElbowPlot(obj_filtered, ndims = 50) + 
  ggtitle("PCA Elbow Plot - Variance Explained")

ggsave(
  filename = "Figures/04_ElbowPlot.png", 
  plot = elbow_plot, 
  width = 8, 
  height = 6, 
  dpi = 300
)

## 3.5 Compute UMAP (Pre-batch correction)
# Run UMAP based on PCA
obj_filtered <- RunUMAP(
  obj_filtered, 
  reduction = "pca", 
  dims = 1:pcs_90,
  seed.use = 1,  
  verbose = FALSE
)

# Visualize batch effect
batch_effect_plot <- DimPlot(
  obj_filtered, 
  reduction = "umap", 
  group.by = "pool", raster = FALSE
) + 
  ggtitle("UMAP: Batch Effect (Pre-Integration)")

ggsave(
  filename = "Figures/05_UMAP_Pre_Integration.png", 
  plot = batch_effect_plot, 
  width = 10, 
  height = 8, 
  dpi = 300
)

### 4. Batch Effect Correction Using Integration Methods
# Evaluation of three integration methods (Harmony, RPCA, and FastMNN) and selection of the best one using LISI metrics.

## 4.1 Harmony Integration

obj_harmony <- obj_filtered

# Run Harmony integration
obj_harmony <- IntegrateLayers(
  object = obj_harmony,
  method = HarmonyIntegration,
  orig.reduction = "pca",
  new.reduction = "harmony",
  verbose = FALSE
)

# Clustering on Harmony reduction
obj_harmony <- FindNeighbors(
  obj_harmony, 
  reduction = "harmony", 
  dims = 1:pcs_90,
  verbose = FALSE
)

obj_harmony <- FindClusters(
  obj_harmony, 
  resolution = 0.3,  
  algorithm = 4,
  random.seed = 1,  
  verbose = FALSE
)

# UMAP on Harmony reduction
obj_harmony <- RunUMAP(
  obj_harmony, 
  reduction = "harmony", 
  dims = 1:pcs_90,
  seed.use = 1,  
  verbose = FALSE
)

# Visualization
harmony_plot <- (
  DimPlot(obj_harmony, reduction = "umap", group.by = "pool") + 
    ggtitle("Batch Effect (Pool)") |
    DimPlot(obj_harmony, reduction = "umap", group.by = "diagnosis") + 
    ggtitle("Diagnosis") |
    DimPlot(obj_harmony, reduction = "umap", label = TRUE) + 
    ggtitle("Clusters")
)

ggsave(
  filename = "Figures/06_Harmony_Integration.png", 
  plot = harmony_plot, 
  width = 18, 
  height = 5, 
  dpi = 300
)

## 4.2 RPCA Integration

obj_rpca <- obj_filtered

obj_rpca <- IntegrateLayers(
  object = obj_rpca, 
  method = RPCAIntegration,
  orig.reduction = "pca", 
  new.reduction = "integrated.rpca",
  verbose = FALSE
)

obj_rpca <- FindNeighbors(
  obj_rpca, 
  reduction = "integrated.rpca", 
  dims = 1:pcs_90,
  verbose = FALSE
)

obj_rpca <- FindClusters(
  obj_rpca, 
  resolution = 0.3, 
  algorithm = 4,
  random.seed = 1,  
  verbose = FALSE
)

obj_rpca <- RunUMAP(
  obj_rpca, 
  reduction = "integrated.rpca", 
  dims = 1:pcs_90,
  seed.use = 1,  
  verbose = FALSE
)

# Visualization
rpca_plot <- (
  DimPlot(obj_rpca, reduction = "umap", group.by = "pool") + 
    ggtitle("Batch Effect (Pool)") |
    DimPlot(obj_rpca, reduction = "umap", group.by = "diagnosis") + 
    ggtitle("Diagnosis") |
    DimPlot(obj_rpca, reduction = "umap", label = TRUE) + 
    ggtitle("Clusters")
)

ggsave(
  filename = "Figures/07_RPCA_Integration.png", 
  plot = rpca_plot, 
  width = 18, 
  height = 5, 
  dpi = 300
)

## 4.3 FastMNN Integration

obj_fastmnn <- obj_filtered

obj_fastmnn <- IntegrateLayers(
  object = obj_fastmnn, 
  method = FastMNNIntegration,
  new.reduction = "integrated.mnn",
  verbose = FALSE
)

obj_fastmnn <- FindNeighbors(
  obj_fastmnn, 
  reduction = "integrated.mnn", 
  dims = 1:pcs_90,
  verbose = FALSE
)

obj_fastmnn <- FindClusters(
  obj_fastmnn, 
  resolution = 0.3, 
  algorithm = 4,
  random.seed = 1,  
  verbose = FALSE
)

obj_fastmnn <- RunUMAP(
  obj_fastmnn, 
  reduction = "integrated.mnn", 
  dims = 1:pcs_90,
  seed.use = 1,  
  verbose = FALSE
)

# Visualization
fastmnn_plot <- (
  DimPlot(obj_fastmnn, reduction = "umap", group.by = "pool") + 
    ggtitle("Batch Effect (Pool)") |
    DimPlot(obj_fastmnn, reduction = "umap", group.by = "diagnosis") + 
    ggtitle("Diagnosis") |
    DimPlot(obj_fastmnn, reduction = "umap", label = TRUE) + 
    ggtitle("Clusters")
)

ggsave(
  filename = "Figures/08_FastMNN_Integration.png", 
  plot = fastmnn_plot, 
  width = 18, 
  height = 5, 
  dpi = 300
)

## 4.4 Evaluate Integration Quality Using LISI Metrics

# LISI (Local Integration Scaled index) is a metric that evaluates:
  
# - Pool LISI: How well batches are mixed (closer to n_batches = better mixing).
# - Cluster LISI: How well cell types are preserved (closer to 1 = better 
# preservation).

# Harmony LISI
lisi_harmony <- compute_lisi(
  X = Embeddings(obj_harmony, "harmony"), 
  meta_data = obj_harmony@meta.data, 
  label_colnames = c("pool", "seurat_clusters")
)

harmony_pool_lisi <- mean(lisi_harmony$pool)
harmony_cluster_lisi <- mean(lisi_harmony$seurat_clusters)

# RPCA LISI
lisi_rpca <- compute_lisi(
  X = Embeddings(obj_rpca, "integrated.rpca"), 
  meta_data = obj_rpca@meta.data, 
  label_colnames = c("pool", "seurat_clusters")
)

rpca_pool_lisi <- mean(lisi_rpca$pool)
rpca_cluster_lisi <- mean(lisi_rpca$seurat_clusters)

# FastMNN LISI
lisi_fastmnn <- compute_lisi(
  X = Embeddings(obj_fastmnn, "integrated.mnn"), 
  meta_data = obj_fastmnn@meta.data, 
  label_colnames = c("pool", "seurat_clusters")
)

fastmnn_pool_lisi <- mean(lisi_fastmnn$pool)
fastmnn_cluster_lisi <- mean(lisi_fastmnn$seurat_clusters)

# Create summary table
lisi_summary <- data.frame(
  Method = c("Harmony", "RPCA", "FastMNN"),
  Pool_LISI = c(harmony_pool_lisi, rpca_pool_lisi, fastmnn_pool_lisi),
  Cluster_LISI = c(harmony_cluster_lisi, rpca_cluster_lisi, fastmnn_cluster_lisi)
)

cat("LISI Metrics Summary\n")
print(lisi_summary)

### 5. Save Integrated Object
# We select the best-integrated object (RPCA in this case) and save it for 
# downstream analysis.

# Select RPCA: best integration method
obj_final <- obj_rpca

# Save the integrated object
saveRDS(obj_final, file = "Results/Integrated_Seurat_Object.rds")


