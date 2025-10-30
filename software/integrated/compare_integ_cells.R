# ==============================================================================
# Compare human-specific cells from various integration methods.
# ==============================================================================

# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(tidyverse)
library(ggalluvial)
library(Seurat)

# set data and output directories
root_dir <- "CBL_scRNAseq"

cca_rds_dir <- file.path("", root_dir, "results/integrated/vladoiu_liam_RL.rds")
harmony_rds_dir <- file.path("", root_dir, "results/integrated/20221003/vladoiu_liam_RL_harmony.rds")
date_dir <- file.path("", root_dir, "results/integrated", format(Sys.Date(), "%Y%m%d"))

if (!dir.exists(date_dir)) {
  dir.create(date_dir)
}

# ------------------------------------------------------------------------------
# import data

cca_srat <- readRDS(cca_rds_dir)
harmony_srat <- readRDS(harmony_rds_dir)

# ------------------------------------------------------------------------------
# reorganize data for plotting

# get metadata from CCA integration
cca_meta <- filter(cca_srat@meta.data, species == "human") %>% 
  select(cell_type = new_cell_type,
         cca_cluster = seurat_clusters) %>% 
  rownames_to_column(var = "cell_id") %>% 
  relocate(cell_id)

# get metadata from harmony integration
harmony_meta <- filter(harmony_srat@meta.data, species == "human") %>% 
  select(cell_type = new_cell_type,
         harmony_cluster = seurat_clusters) %>% 
  rownames_to_column(var = "cell_id") %>% 
  relocate(cell_id)

# join both and convert to alluvial format
meta_df <- full_join(x = cca_meta,
                     y = harmony_meta,
                     by = c("cell_id", "cell_type")) %>% 
  count(cca_cluster, harmony_cluster, cell_type, name = "frequency")

# ------------------------------------------------------------------------------
# make alluvial plot

plot_alluvial <- function(data) {
  plt <- ggplot(data,
                aes(y = frequency, axis1 = cca_cluster, axis2 = harmony_cluster)) + 
    geom_alluvium(aes(fill = cell_type), width = 0.1) + 
    geom_stratum(width = 0.1) + 
    geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 3) + 
    labs(x = "clustering", y = "number of cells") + 
    scale_x_discrete(limits = c("CCA cluster", "harmony cluster"),
                     expand = expansion(mult = 0.1))
  
  return(plt)
}

ggsave(
  "alluvial_allcells.png",
  plot = plot_alluvial(meta_df),
  path = date_dir,
  width = 5,
  height = 7.5,
  units = "in",
  dpi = 600
)

ggsave(
  "alluvial_UBCs.png",
  plot = plot_alluvial(meta_df[meta_df$cca_cluster %in% c(7, 19, 20) |
                                 meta_df$harmony_cluster %in% c(5, 16),]),
  path = date_dir,
  width = 5,
  height = 5,
  units = "in",
  dpi = 600
)
