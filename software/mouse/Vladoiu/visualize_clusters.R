# ==============================================================================
# Takes a Seurat object, finds the clusters, and plots t-SNE/UMAP
# ==============================================================================

# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
# setwd("/CBL_scRNAseq/software/")

library(tidyverse)
library(Seurat)

# set data and output directories
root_dir <- "CBL_scRNAseq"
data_dir <- file.path("/isilon", root_dir, "data/mouse/Vladoiu")
out_dir <- file.path("", root_dir, "results/mouse/Vladoiu")
date_dir <- file.path(out_dir, format(Sys.Date(),"%Y%m%d"))

if (!dir.exists(date_dir)) {
  dir.create(date_dir)
}

# load Seurat object from RDS file
srat <- readRDS(file.path(out_dir, "merged_seurat.rds"))

# get variable features
var_feat <- VariableFeatures(srat)

# run PCA, check PCs
srat <- RunPCA(srat, features = var_feat, npcs = 100)
DimHeatmap(srat, dims = 1:15, cells = 500)
ElbowPlot(srat, ndims = 100)

# set ndims for clustering and dimensional reductions (based on elbow plot)
ndims <- 30

# cluster the cells
srat <- FindNeighbors(srat, dims = 1:ndims)
srat <- FindClusters(srat)

# run t-SNE
srat <- RunTSNE(srat, dims = 1:ndims, num_threads = 32)
DimPlot(srat, reduction = "tsne", label = TRUE, repel = TRUE) + NoLegend()
ggsave("tsne_seurat_clusters.png", path = date_dir, width = 4, height = 4, units = "in", dpi = 600)
DimPlot(srat, reduction = "tsne", group.by = "cell_type", label = TRUE, repel = TRUE) + NoLegend()
ggsave("tsne_cell_type_annot.png", path = date_dir, width = 4, height = 4, units = "in", dpi = 600)

# run UMAP
srat <- RunUMAP(srat, dims = 1:ndims)
DimPlot(srat, reduction = "umap", label = TRUE, repel = TRUE) + NoLegend()
ggsave("umap_seurat_clusters.png", path = date_dir, width = 4, height = 4, units = "in", dpi = 600)
DimPlot(srat, reduction = "umap", group.by = "cell_type", label = TRUE, repel = TRUE) + NoLegend()
ggsave("umap_cell_type_annots.png", path = date_dir, width = 4, height = 4, units = "in", dpi = 600)

# save the clustered Seurat object
srat_file <- file.path(out_dir, "merged_seurat.rds")
message(sprintf("Saving the clustered object as %s", srat_file))
saveRDS(srat, file = srat_file)

# ------------------------------------------------------------------------------
# rename Vladoiu annotations and make plot

source(file.path("", root_dir, "software/mouse/Vladoiu/add_annotations.R"))
srat <- label_vladoiu_cells(srat)

annots <- table(srat$cell_type, srat$seurat_clusters) %>% 
  as.data.frame(.)
colnames(annots) <- c("annot", "cluster", "freq")

# plot number of annotated mouse cells in each cluster
ggplot(annots, aes(x = cluster, y = freq, fill = annot)) + 
  geom_col() + 
  labs(y = "number") + 
  theme_light() + 
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) + 
  # use dittoSeq colour palette: 
  # https://github.com/dtm2451/dittoSeq/blob/master/R/dittoColors.R
  scale_fill_manual(values = c("#E69F00", "#56B4E9", "#009E73", "#F0E442",
                               "#0072B2", "#D55E00", "#CC79A7", "#666666",
                               "#AD7700", "#1C91D4", "#007756", "#D5C711",
                               "#005685", "#A04700", "#B14380", "#4D4D4D",
                               "#FFBE2D", "#80C7EF", "#00F6B3", "#F4EB71",
                               "#06A5FF", "#FF8320", "#D99BBD", "#8C8C8C",
                               "#FFCB57", "#9AD2F2", "#2CFFC6", "#F6EF8E",
                               "#38B7FF", "#FF9B4D", "#E0AFCA", "#A3A3A3",
                               "#8A5F00", "#1674A9", "#005F45", "#AA9F0D",
                               "#00446B", "#803800", "#8D3666", "#3D3D3D"))
ggsave("mouse_cells_per_cluster.pdf", path = date_dir, 
       width = 20, height = 8, units = "in")
