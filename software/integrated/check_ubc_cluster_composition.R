# ==============================================================================
# Check composition of the clusters (e.g., by age, by cell type, etc.)
# ==============================================================================

# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(tidyverse)
library(Seurat)

# set data and output directories
root_dir <- "CBL_scRNAseq"
out_dir <- file.path(paste0("/", root_dir), "results/integrated")
date_dir <- file.path(out_dir, "20220726")

# ------------------------------------------------------------------------------
# load data

srat <- readRDS(file.path("", root_dir, "results/human/Aldinger/UBC_seurat.rds"))

# ------------------------------------------------------------------------------
# visualize clusters as bar graphs

source(file.path("", root_dir, "software/utilities/cluster_barplot.R"))

all_split <- c("integ_clusters", rep(c("age", "new_cell_type"), 2))
all_group <- c(rep("seurat_clusters", 3), rep("integ_clusters", 2))

all_plt <- apply(
  X = cbind(all_split, all_group),
  MARGIN = 1,
  FUN = function(X) {
    plt <- cluster_barplot(
      srat,
      split.by = X[1],
      group.by = X[2],
      position = "stack"
    )
    
    ggsave(
      filename = paste0(X[1], "-by-", X[2], ".png"),
      plot = plt,
      path = date_dir,
      width = 5,
      height = 3,
      units = "in",
      dpi = 300
    )
    
    return(plt)
  }
) %>% 
  patchwork::wrap_plots(., ncol = 3)

ggsave(
  filename = "combined_plots.png",
  plot = all_plt,
  path = date_dir,
  width = 12,
  height = 6,
  units = "in",
  dpi = 300
)
