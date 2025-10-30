# ==============================================================================
# Visualize and explore integrated data
# ==============================================================================

# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(argparse)
library(tidyverse)
library(Seurat)

# parse command line arguments
parser <- ArgumentParser()
parser$add_argument("--out_dir", default = NULL)
parser$add_argument("--date_dir", default = NULL)
parser$add_argument("--srat_rds", default = NULL)

args <- parser$parse_args()

# set data and output directories
root_dir <- "CBL_scRNAseq"
human_data_dir <- file.path("/isilon", root_dir, "data/human")
mouse_data_dir <- file.path("/isilon", root_dir, "data/mouse")

if (is.null(args$out_dir)) {
  out_dir <- file.path(paste0("/", root_dir), "results/integrated")
} else {
  out_dir <- args$out_dir
}

if (is.null(args$date_dir)) {
  date_dir <- file.path(out_dir, format(Sys.Date(), "%Y%m%d"))
} else {
  date_dir <- args$date_dir
}

if (is.null(args$srat_rds)) {
  srat_rds <- file.path(out_dir, "vladoiu_liam_RL.rds")
} else {
  srat_rds <- args$srat_rds
}

message(sprintf("Reading Seurat object from: %s", srat_rds))
message(sprintf("Saving results to: %s", date_dir))

if (!dir.exists(date_dir)) {
  dir.create(date_dir)
}

# ==============================================================================
# import data
# ==============================================================================

srat <- readRDS(srat_rds)

# ==============================================================================
# run basic visual checks
# ==============================================================================

source(file.path("", root_dir, "software/utilities/initial_visual_check.R"))

# ==============================================================================
# based on Liam's annotations, what human cells are in each cluster?
# ==============================================================================

human_annot <- table(srat$new_cell_type, srat$seurat_clusters) %>% 
  as.data.frame(.)
colnames(human_annot) <- c("annot", "cluster", "freq")

# plot number of annotated human cells in each cluster
ggplot(human_annot, aes(x = cluster, y = freq, fill = annot)) + 
  geom_col() + 
  labs(y = "number") + 
  theme_light() + 
  scale_y_continuous(expand = expansion(mult = c(0, 0.05)))
ggsave("human_cells_per_cluster.pdf", path = date_dir, 
       width = 10, height = 5, units = "in")

# ==============================================================================
# based on annotations from Vladoiu paper, what mouse cells are in each cluster?
# ==============================================================================

# change mouse cell types to factor and rename
source(file.path("", root_dir, "software/mouse/Vladoiu/add_annotations.R"))
srat <- label_vladoiu_cells(srat)

mouse_annot <- table(srat$mouse_cell_type, srat$seurat_clusters) %>% 
  as.data.frame(.)
colnames(mouse_annot) <- c("annot", "cluster", "freq")

# plot number of annotated mouse cells in each cluster
ggplot(mouse_annot, aes(x = cluster, y = freq, fill = annot)) +
  geom_col() +
  labs(y = "number") +
  theme_light() +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  scale_fill_manual(
    values = c(
      # RL/UBC lineage
      scales::brewer_pal(palette = "Reds")(4),
      # granule cell (progenitors)
      scales::brewer_pal(palette = "Blues")(3),
      # other glutamatergic
      scales::brewer_pal(palette = "YlOrBr")(3),
      # stem cells
      scales::brewer_pal(palette = "BuPu")(5),
      # GABAergic lineage
      scales::brewer_pal(palette = "Greens")(7),
      # glial/non-neuronal cells
      scales::brewer_pal(palette = "PuRd")(9),
      # NA
      "#DFDFDF"
    )
  )
ggsave("mouse_cells_per_cluster.pdf", path = date_dir, 
       width = 15, height = 5, units = "in")

# ==============================================================================
# which timepoints do each of the cells/clusters come from?
# ==============================================================================

cell_ages <- table(srat$orig.ident, srat$seurat_clusters) %>% 
  as.data.frame(.)
colnames(cell_ages) <- c("timepoint", "cluster", "freq")

# plot cells from each timepoint
ggplot(cell_ages, aes(x = cluster, y = freq, fill = timepoint)) + 
  geom_col() + 
  labs(y = "number") + 
  theme_light() + 
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) + 
  scale_fill_brewer(palette = "Set3")
ggsave("timepoints_per_cluster_bar.pdf", path = date_dir,
       width = 10, height = 5, units = "in")

DimPlot(srat, group.by = "orig.ident") + 
  scale_colour_brewer(palette = "Set3")
ggsave("timepoints_per_cluster_umap.pdf", path = date_dir,
       width = 10, height = 7.5, units = "in")
