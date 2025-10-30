# ==============================================================================
# Isolate RL-derived cells from the merged Seurat object
# ==============================================================================

# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(tidyverse)
library(Seurat)

# set data and output directories
root_dir <- "CBL_scRNAseq"
data_dir <- file.path("/isilon", root_dir, "data/mouse/Vladoiu")
out_dir <- file.path("", root_dir, "results/mouse/Vladoiu")
date_dir <- file.path(out_dir, format(Sys.Date(), "%Y%m%d"))

if (!dir.exists(date_dir)) {
  dir.create(date_dir)
}

# load data
full_srat <- readRDS(file.path(out_dir, "merged_seurat.rds"))

# isolate RL-derived cells
# info in the lab notebook (2022-04-21)
# updated 2022-09-11
rlip_clusters <- c(0:1, 3, 6, 9, 11, 13, 16, 21:22, 34)
rlip_srat <- subset(full_srat, idents = rlip_clusters)

# show RL-derived cells on UMAP of all cells
rlip_cells <- WhichCells(full_srat, idents = rlip_clusters)
viz_rlip <- DimPlot(full_srat, cells.highlight = rlip_cells, 
                    cols.highlight = "orange", sizes.highlight = 0.01) + 
  NoLegend()
viz_clusters <- DimPlot(full_srat, label = TRUE, repel = TRUE) + 
  NoLegend()
viz_rlip + viz_clusters
ggsave("glutamatergic_clusters.pdf", 
       path = date_dir, 
       width = 8, height = 4, units = "in")

# save subsetted RL-derived cells as new RDS object
saveRDS(rlip_srat, file = file.path(out_dir, "merged_seurat_RLonly.rds"))
