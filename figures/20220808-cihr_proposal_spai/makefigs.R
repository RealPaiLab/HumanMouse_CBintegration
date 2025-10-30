# ------------------------------------------------------------------------------
# Figures requestd by Shraddha on 2022-08-08 for CIHR proposal
# ------------------------------------------------------------------------------

# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(Seurat)
library(tidyverse)
library(patchwork)

# function for showing the composition of each cluster
source("/CBL_scRNAseq/software/utilities/cluster_barplot.R")

# import Seurat object
srat <- readRDS("/CBL_scRNAseq/results/integrated/vladoiu_liam_RL.rds")

# colour of human vs mouse
species_cols <- c("#007FA3", "grey")

# UMAPs
umap_cell <- DimPlot(srat, group.by = "seurat_clusters", label = TRUE, repel = TRUE) + 
  NoLegend()
umap_species <- DimPlot(srat, label = FALSE, group.by = "species", pt.size = 0.01, 
                        cells.highlight = rownames(srat@meta.data)[srat$species == "human"],
                        sizes.highlight = 0.01) + 
  scale_colour_manual(labels = c("mouse", "human"),
                      values = rev(species_cols))
umaps <- umap_cell + umap_species + plot_annotation(tag_levels = "A") & labs(title = NULL)
ggsave("umaps.png", plot = umaps, width = 8, height = 3, units = "in")

# UMAP by human cell type
srat$new_cell_type <- case_when(
  is.na(srat$new_cell_type) == TRUE ~ "NA",
  TRUE ~ srat$new_cell_type
)
umap_human_cells <- DimPlot(srat, group.by = "new_cell_type", 
                            order = rev(names(table(srat$new_cell_type))[-6])) + 
  labs(title = NULL) + 
  scale_colour_manual(
    values = c("grey", hcl(h = seq(15, 375, length = 8), c = 100, l = 65)[1:7])
  )
ggsave("umap_human_cells.png", width = 5, height = 3, units = "in")

# barplot
bar_species <- cluster_barplot(srat, split.by = "species") + 
  labs(x = "cell type clusters") + 
  scale_fill_manual(values = species_cols)
ggsave("barplot_species.png", plot = bar_species, width = 6, height = 3, units = "in")
