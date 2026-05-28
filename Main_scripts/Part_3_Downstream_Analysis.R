# ==============================================================================
# Script: Part_3_Downstream_Analysis.R
# Project: Single-cell transcriptomic profiling of peripheral immune cells in 
#          genetically stratified Parkinson’s disease and progressive supranuclear 
#          palsy
# Description: This script performs cell type abundance analysis, pseudobulk 
#              differential expression (DEG) analysis using DESeq2, Venn diagram
#              overlap visualization, and hierarchical clustering heatmaps 
#              to reveal lineage-specific gene expression patterns.
# ==============================================================================

# Create output directories
dir.create("Results/Venn_HC", showWarnings = FALSE)
dir.create("Figures/Venn_Pathology", showWarnings = FALSE)
dir.create("Figures/Heatmaps", showWarnings = FALSE)

# Load Required Libraries
# sc-RNAseq analysis
library(Seurat)

# Differential Expression Analysis
library(DESeq2)     

# Visualization
library(pheatmap)        
library(ggvenn)         
library(ggpubr)  
library(ggplot2)
library(patchwork)
library(RColorBrewer)

# Statistical testing
library(rstatix)    

# Data handling
library(openxlsx)  
library(tidyverse)
library(dplyr)

### 1. Load Annotated Object from Part 2.
obj_final <- readRDS("Results/Final_Annotated_Seurat_Object.rds")

# Define cell populations for analysis
tipos_principales <- c(
  "T_CD4_cells", "T_CD8_cells", "NK_cells", "B_cells", 
  "Monocytes", "Dendritic_cells"
)

# Cell types used in pseudobulk aggregation
tipos_pseudobulk <- c(
  "Monocytes", "Dendritic-cells", 
  "T-CD4-cells", "T-CD8-cells", "NK-cells","B-cells"
)

### 1. Relative Abundance Analysis
# Analysis of the proportions of each cell type across diagnostic groups.

## 1.1 Calculate Cell Type Proportions
# Calculate counts per cell type per sample
datos_por_muestra <- obj_final@meta.data %>%
  filter(!is.na(final_annotation), final_annotation %in% tipos_principales) %>%
  group_by(id_tecnico, comparison_group, final_annotation) %>%
  summarise(n = n(), .groups = "drop") %>%
  # Calculate percentage within each sample
  group_by(id_tecnico) %>%
  mutate(porcentaje = n / sum(n) * 100) %>%
  ungroup() %>%
  mutate(final_annotation = factor(final_annotation, levels = tipos_principales))

## 1.2 Visualize Mean Proportions
# Calculate mean proportions per diagnosis for visualization
datos_grafico <- datos_por_muestra %>%
  group_by(comparison_group, final_annotation) %>%
  summarise(porcentaje_medio = mean(porcentaje), .groups = "drop")

# Create stacked bar plot
abundance_plot <- ggplot(
  datos_grafico, 
  aes(x = comparison_group, y = porcentaje_medio, fill = final_annotation)
) +
  geom_col(width = 0.6, color = "white", size = 0.3) +
  scale_fill_manual(
    values = c(
      "T_CD4_cells" = "#6a99b5",
      "T_CD8_cells" = "#7fa998",
      "NK_cells" = "#f3a661",
      "B_cells" = "#b39bc8",
      "Monocytes" = "#df7a76",
      "Dendritic_cells" = "#f2cc8f"
    )
  ) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 100.1)) +
  theme_bw(base_size = 14) + 
  theme(plot.title = element_text(size = 25, face = "bold", hjust = 0.5, margin = margin(b = 20)),
        axis.title.y = element_text(size = 20, face = "bold", margin = margin(r = 10)),
        axis.title.x = element_text(size = 20, face = "bold", margin = margin(r = 10)),
        axis.text.x = element_text(color = "black", size = 20, face = "bold", angle = 45, hjust = 1),
        axis.text.y = element_text(color = "black", size = 20),
        
        # Leyenda: Más grande para que el texto sea claro
        legend.title = element_text(size = 18, face = "bold"),
        legend.text = element_text(size = 18),
        legend.position = "right",
        plot.margin = margin(10, 10, 10, 10)
  ) + labs(
    title = "Relative Cell Type Abundance Across Diagnostic Groups",
    y = "Mean Relative Proportion (%)", x = "Diagnostic groups",
    fill = "Cell Type"
  )

ggsave(
  "Figures/14_RelativeAbundance_BarPlot.png",
  plot = abundance_plot,
  width = 20, 
  height = 8, 
  dpi = 300
)

## 1.3 Statistical Testing (Wilcoxon)
# Test if cell type proportions differ significantly between diagnostic groups.
resultados_wilcox <- list()

for (celltype in tipos_principales) {
  # Get data for this cell type
  datos_celltype <- datos_por_muestra %>% 
    filter(final_annotation == celltype)
  
  # Perform pairwise Wilcoxon test with multiple testing correction
  test_result <- pairwise.wilcox.test(
    x = datos_celltype$porcentaje, 
    g = datos_celltype$comparison_group, 
    p.adjust.method = "BH"
  )
  
  # Convert to data frame format
  tabla_resultado <- as.data.frame(test_result$p.value) %>%
    rownames_to_column(var = "Group_1") %>%
    pivot_longer(
      cols = -Group_1, 
      names_to = "Group_2", 
      values_to = "adjusted p-value"
    ) %>%
    filter(!is.na(p_adj)) %>%
    mutate(
      Cell_Type = celltype,
      p_adj = round(p_adj, 4),
      Significant = ifelse(p_adj < 0.05, "Yes", "No")
    )
  
  resultados_wilcox[[celltype]] <- tabla_resultado
}

# Combine all results
wilcox_abundance <- bind_rows(resultados_wilcox) %>%
  select(Cell_Type, Group_1, Group_2, p_adj, Significant)

# Save results
write.xlsx(
  wilcox_abundance,
  "Results/Relative_Abundance_Wilcoxon_Results.xlsx",
  overwrite = TRUE
)

# 2. Pseudobulk Differential Expression Analysis
# Pseudobulk aggregation combines cells from the same cell type in the
# same sample into a single "bulk" sample, which is then tested for
# differential expression using DESeq2.

## 2.1 Create Pseudobulk Matrices
# Aggregate cells by cell type, diagnosis, and sample ID
pseudo_counts <- AggregateExpression(
  obj_final,
  group.by = c("final_annotation", "comparison_group", "id_tecnico"),
  assays = "RNA",
  return.seurat = TRUE
)

# Create grouping identifier
pseudo_counts$celltype_group <- paste(
  pseudo_counts$final_annotation,
  pseudo_counts$comparison_group,
  sep = "_"
)

Idents(pseudo_counts) <- "celltype_group"

## 2.2 Differential Expression Analysis
# Define comparisons to test
diagnosis <- c("sporadic-PD","Genetic-PD", "PSP", "HC")
comparisons <- combn(diagnosis, 2, simplify = FALSE)

# Store results
de_results <- list()

for (celltype in tipos_pseudobulk) {
  for (comp in comparisons) {
    group1 <- comp[1]
    group2 <- comp[2]
    
    # Create identifiers for this comparison
    id1 <- paste(celltype, group1, sep = "_")
    id2 <- paste(celltype, group2, sep = "_")
    
    # Check if both groups exist in data
    if (!all(c(id1, id2) %in% Idents(pseudo_counts))) {
      next
    }
    
    # Need at least 3 samples per group
    n1 <- sum(Idents(pseudo_counts) == id1)
    n2 <- sum(Idents(pseudo_counts) == id2)
    
    if (n1 < 3 | n2 < 3) {
      next
    }
    
    cat(sprintf("  Testing %s: %s (n=%d) vs %s (n=%d)\n",
                celltype, group1, n1, group2, n2))
    
    # Perform differential expression test
    degs <- tryCatch({
      FindMarkers(
        pseudo_counts,
        ident.1 = id1,
        ident.2 = id2,
        test.use = "DESeq2",
        verbose = FALSE
      )
    }, error = function(e) {
      message(sprintf("    Error: %s", e$message))
      return(NULL)
    })
    
    # Extract and filter significant genes
    if (!is.null(degs) && nrow(degs) > 0) {
      sig_degs <- degs %>%
        rownames_to_column(var = "gene") %>%
        # Filter by statistical significance and effect size
        filter(p_val_adj < 0.05, abs(avg_log2FC) > 0.5) %>%
        mutate(
          comparison = paste0(group1, "_vs_", group2),
          cell_type = celltype,
          group_1 = group1,
          group_2 = group2
        ) %>%
        arrange(p_val_adj)
      
      if (nrow(sig_degs) > 0) {
        de_results[[paste0(celltype, "_", group1, "_", group2)]] <- sig_degs
      }
    }
  }
}

# Combine all DE results
if (length(de_results) > 0) {
  final_de_table <- bind_rows(de_results)
  
  write.xlsx(
    final_de_table,
    "Results/Pseudobulk_DEG_Results_AllComparisons.xlsx",
    overwrite = TRUE
  )
  
  cat(sprintf("  Total significant DEGs found: %d\n", nrow(final_de_table)))
} 

### 3. Venn Diagrams - DEGs Overlap Analysis
# Venn diagrams show overlap in differentially expressed genes between different 
# comparisons.

## 3.1 Venn Diagrams: Pathologies vs Healthy Control
if (nrow(final_de_table) > 0) {
  cat("Creating Venn diagrams for pathologies vs HC\n")
  
  lineages_config <- list(
    Myeloid = list(
      celltypes = c("Monocytes", "Dendritic-cells"),
      colors = c("#FDF2E9", "#E67E22", "#D35400")
    ),
    Lymphoid = list(
      celltypes = c("T-CD4-cells", "T-CD8-cells", "NK-cells", "B-cells"),
      colors = c("#F4ECF7", "#8E44AD", "#5B2C6F")
    )
  )
  
  pathology_comparisons <- c("sporadic-PD_vs_HC", "Genetic-PD_vs_HC", "PSP_vs_HC")
  
  for (lineage_name in names(lineages_config)) {
    lineage_info <- lineages_config[[lineage_name]]
    celltypes <- lineage_info$celltypes
    colors <- lineage_info$colors
    
    # Extract genes for each pathology
    genes_by_pathology <- list()
    
    for (comp in pathology_comparisons) {
      pathology_name <- gsub("_vs_HC", "", comp)
      
      genes_union <- final_de_table %>%
        filter(cell_type %in% celltypes, comparison == comp) %>%
        pull(gene) %>%
        unique()
      
      if (length(genes_union) > 0) {
        genes_by_pathology[[pathology_name]] <- genes_union
      }
    }
    
    # Create Venn diagram if we have at least 2 pathologies
    if (length(genes_by_pathology) >= 2) {
      venn_plot <- ggvenn(
        genes_by_pathology,
        columns = names(genes_by_pathology),
        fill_color = colors[1:length(genes_by_pathology)],
        stroke_color = "grey40",
        stroke_size = 0.8,
        set_name_size = 5,
        text_size = 5,
        show_percentage = TRUE
      ) +
        labs(
          title = paste("DEG Overlap:", lineage_name, "Lineage")
        ) +
        theme_void(base_family = "Helvetica") +
        theme(
          plot.title = element_text(face = "bold", size = 12, hjust = 0.5, margin = margin(b = 5)),
          plot.margin = margin(15, 15, 15, 15)
        )
      
      ggsave(
        filename = paste0("Figures/Venn_HC/Venn_", lineage_name, "_vs_HC.png"),
        plot = venn_plot,
        width = 6,
        height = 5,
        dpi = 300,
        bg = "white"
      )
    }
  }
}

## 3.2 Venn Diagrams: Cross-diagnostic comparisons
if (nrow(final_de_table) > 0) {
  cat("Creating Venn diagrams for comparisons between pathologies...\n")
  
  pathology_comparisons <- c(
    "sporadic-PD_vs_Genetic-PD",
    "sporadic-PD_vs_PSP",
    "Genetic-PD_vs_PSP"
  )
  
  for (lineage_name in names(lineages_config)) {
    lineage_info <- lineages_config[[lineage_name]]
    celltypes <- lineage_info$celltypes
    colors <- lineage_info$colors
    
    # Extract genes for each comparison
    genes_by_comparison <- list()
    
    for (comp in pathology_comparisons) {
      genes_union <- final_de_table %>%
        filter(cell_type %in% celltypes, comparison == comp) %>%
        pull(gene) %>%
        unique()
      
      if (length(genes_union) > 0) {
        genes_by_comparison[[comp]] <- genes_union
      }
    }
    
    # Create Venn diagram if we have at least 2 comparisons
    if (length(genes_by_comparison) >= 2) {
      formatted_names <- names(genes_by_comparison) %>%
        gsub("_", " ", .) %>% gsub(" vs ", " vs\n", .)
      
      names(genes_by_comparison) <- formatted_names
      
      venn_plot <- ggvenn(
        genes_by_comparison,
        columns = names(genes_by_comparison),
        fill_color = colors[1:length(genes_by_comparison)],
        stroke_color = "grey40",
        stroke_size = 0.8,
        set_name_size = 6,    
        text_size = 7,
        show_percentage = TRUE
      ) + 
        theme_void() +
        theme(plot.margin = margin(t = 30, r = 20, b = 20, l = 20))
      
      title_plot <- ggplot() + 
        annotate("text", x = 1, y = 1, label = paste("DEGs Overlap:", lineage_name, "Lineage"), 
                 fontface = "bold", size = 7) + 
        theme_void()
      final_plot <- title_plot / venn_plot + plot_layout(heights = c(0.1, 1))
      
      ggsave(
        filename = paste0("Figures/Venn_Pathology/Venn_", lineage_name, "_PathologyComparisons.png"),
        plot = final_plot,
        width = 8,
        height = 8,
        dpi = 300,
        bg = "white"
      )
    }
  }
}

### 4. Heatmaps of DEGs: hierarchical clustering by samples
# Create heatmaps showing expression of differentially expressed genes.

## 4.1 Prepare Data for Heatmaps
# Load DEGs table
final_de_table <- read.xlsx("Results/Pseudobulk_DEG_Results_AllComparisons.xlsx")

# Diagnosis color palette
annotation_colors <- list(
  Diagnosis = c(
    "HC"           = "#4CAF50",
    "sporadic-PD"  = "#E64B35",
    "Genetic-PD" = "#4DBBD5",
    "PSP"          = "#9B59B6"
  )
)

# Heatmap color scale (blue = low, white = neutral, red = high)
heatmap_colors <- colorRampPalette(
  c("#2166AC", "#92C5DE", "#F7F7F7", "#F4A582", "#D6604D")
)(100)

# Manual mapping of pseudobulk sample IDs to diagnosis groups
# (pseudobulk uses '-' separator from AggregateExpression)
sample_diagnosis <- c(
  "HC-A"    = "HC",           "HC-B"    = "HC",
  "HC-C"    = "HC",           "HC-D"    = "HC",
  "sPD-A"   = "sporadic-PD",  "sPD-B"   = "sporadic-PD",
  "sPD-C"   = "sporadic-PD",  "sPD-D"   = "sporadic-PD",
  "GBA-B"   = "Genetic-PD", "GBA-D"   = "Genetic-PD",
  "LRRK2-A" = "Genetic-PD", "LRRK2-C" = "Genetic-PD",
  "PSP-A"   = "PSP",          "PSP-B"   = "PSP",
  "PSP-C"   = "PSP",          "PSP-D"   = "PSP"
)

# Lineage definitions
lineages <- list(
  Myeloid = list(
    celltypes_seurat = c("Monocytes", "Dendritic_cells"),
    celltypes_table  = c("Monocytes", "Dendritic-cells"),
    outfile          = "Figures/Heatmaps/Heatmap_Myeloid.png"
  ),
  Lymphoid = list(
    celltypes_seurat = c("T_CD4_cells", "T_CD8_cells", "NK_cells", "B_cells"),
    celltypes_table  = c("T-CD4-cells", "T-CD8-cells", "NK-cells", "B-cells"),
    outfile          = "Figures/Heatmaps/Heatmap_Lymphoid.png"
  )
)


## 4.2 Heatmap function
make_heatmap <- function(obj, df, celltypes_seurat, celltypes_table,
                         lineage_name, outfile) {
  
  # Step 1: Extract DEG gene list 
  genes_sig <- df %>%
    filter(cell_type %in% celltypes_table) %>%
    pull(gene) %>%
    unique()
  
  message(sprintf("[%s] %d unique DEGs extracted.", lineage_name, length(genes_sig)))
  if (length(genes_sig) < 2) {
    message("fewer than 2 genes.")
    return(invisible(NULL))
  }
  
  # Step 2: Pseudobulk per sample 
  Idents(obj) <- "final_annotation"
  sub <- subset(obj, idents = celltypes_seurat)
  
  pb <- AggregateExpression(
    sub,
    group.by      = "id_tecnico",
    assays        = "RNA",
    return.seurat = TRUE
  )
  
  # Step 3: Normalize 
  pb <- NormalizeData(pb,
                      normalization.method = "LogNormalize",
                      scale.factor = 10000,
                      verbose = FALSE)
  
  exp_mat <- GetAssayData(pb, assay = "RNA", layer = "data")
  
  # Step 4: Filter to DEG genes 
  genes_present <- intersect(genes_sig, rownames(exp_mat))
  message(sprintf("[%s] %d / %d DEGs found in expression matrix.",
                  lineage_name, length(genes_present), length(genes_sig)))
  
  if (length(genes_present) < 2) {
    message("fewer than 2 DEGs found in matrix.")
    return(invisible(NULL))
  }
  
  exp_mat <- as.matrix(exp_mat[genes_present, ])
  
  # Remove genes with >50% zeros
  pct_zeros <- rowMeans(exp_mat == 0)
  exp_mat   <- exp_mat[pct_zeros <= 0.50, ]
  n_genes <- nrow(exp_mat)
  
  if (n_genes < 2) {
    message("fewer than 2 genes after filtering.")
    return(invisible(NULL))
  }
  
  z_mat <- t(scale(t(exp_mat)))
  z_mat[is.nan(z_mat)] <- 0
  z_mat[is.na(z_mat)]  <- 0
  z_mat[z_mat >  2]    <-  2  
  z_mat[z_mat < -2]    <- -2  
  
  muestras    <- colnames(z_mat)
  sample_meta <- data.frame(
    Diagnosis = factor(sample_diagnosis[muestras],
                       levels = c("HC", "sporadic-PD", "Genetic-PD", "PSP")),
    row.names = muestras
  )
  
  # Bold column labels
  bold_colnames <- lapply(muestras, function(x) bquote(bold(.(x))))
  
  # Plot
  pheatmap(
    mat                      = z_mat,
    color                    = heatmap_colors,
    scale                    = "none",
    breaks                   = seq(-2, 2, length.out = 101),
    
    # Clustering
    cluster_cols             = TRUE,
    clustering_distance_cols = "correlation",
    clustering_method        = "ward.D2",
    cluster_rows             = TRUE,
    clustering_distance_rows = "correlation",
    
    # Annotations
    annotation_col           = sample_meta,
    annotation_colors        = annotation_colors,
    
    # Labels
    show_rownames            = (n_genes <= 80),
    show_colnames            = TRUE,
    labels_col               = as.expression(bold_colnames),
    fontsize                 = 18,
    fontsize_col             = 15,
    fontsize_row             = 15,
    angle_col                = "45",
    
    # Aesthetics
    border_color             = "#EFEFEF",
    main                     = sprintf("%s Lineage",
                                       lineage_name, n_genes),
    # Save
    filename                 = outfile,
    width                    = 10,
    height                   = 12
  )
}

## 4.3 Generate Hierarchical Clustering Heatmaps
for (lin in names(lineages)) {
  message(sprintf("\n===== %s LINEAGE =====", toupper(lin)))
  make_heatmap(
    obj              = obj_final,
    df               = final_de_table,
    celltypes_seurat = lineages[[lin]]$celltypes_seurat,
    celltypes_table  = lineages[[lin]]$celltypes_table,
    lineage_name     = lin,
    outfile          = lineages[[lin]]$outfile
  )
}
