# ==============================================================================
# Supplementary figure of manuscript showing UMAPs of the full cerebellum after
# merging the human and mouse datasets (without integration). This script should
# be run using the scrnaseq_env conda environment.
# ==============================================================================

library(tidyverse)
library(patchwork)
library(Seurat)

source("./utils.R")

out_dir <- file.path(
  "/.mounts/labs/pailab/private/projects/HumanMouseUBC/figures/SuppFig_merged_unintegrated",
  format(Sys.Date(), "%Y%m%d")
)
if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

# colour palette
my_pals <- get_custom_pals()

# ------------------------------------------------------------------------------
# UMAPs of the merged datasets

# load merged Seurat object
srat_qs <- get_srat_paths()
srat <- load_srat(srat_qs["merge"]) %>% pluck(1)

srat$dataset_name <- str_remove(
  string = srat$dataset_name,
  pattern = "full_cerebellum_"
)

# plot UMAP by species
umap_species <- DimPlot(
  srat,
  reduction = "umap",
  cells.highlight = WhichCells(srat, expression = species == "human"),
  sizes.highlight = 0.01,
  label = FALSE,
  raster = FALSE
) + 
  labs(title = "Merged human+mouse (species)", x = "UMAP 1", y = "UMAP 2") + 
  scale_colour_manual(
    labels = c("mouse", "human"),
    values = my_pals$species
  ) + 
  theme_classic2() + 
  theme(
    legend.position = "bottom",
    legend.justification = c(0.5, 1)
  )
ggsave(
  filename = "merge_full_species.pdf",
  plot = umap_species,
  path = out_dir,
  width = 5,
  height = 5.5,
  units = "in"
)

# plot UMAP by dataset
umap_dataset <- DimPlot(
  srat,
  reduction = "umap",
  group.by = "dataset_name",
  label = FALSE,
  raster = FALSE
) + 
  labs(title = "Merged human+mouse (dataset)", x = "UMAP 1", y = "UMAP 2") + 
  scale_colour_manual(
    labels = \(x) {str_replace(string = x, pattern = "_(.*)", replacement = " (\\1)")},
    values = alpha(scales::pal_hue()(4), alpha = 0.2)
  ) + 
  guides(colour = guide_legend(
    override.aes = list(size = 3, alpha = 1),
    nrow = 2
  )) + 
  theme_classic2() + 
  theme(
    legend.position = "bottom",
    legend.justification = c(0.5, 1)
  )
ggsave(
  filename = "merge_full_dataset.pdf",
  plot = umap_dataset,
  path = out_dir,
  width = 5,
  height = 5.5,
  units = "in"
)

# ------------------------------------------------------------------------------
# final figure

.plt <- wrap_plots(
  umap_species,
  umap_dataset
) + 
  plot_layout(ncol = 2)

walk(
  .x = c("png", "pdf"),
  .f = \(device) {
    ggsave(
      filename = paste("combined", device, sep = "."),
      plot = .plt,
      path = out_dir,
      width = 7.5,
      height = 4.5,
      units = "in",
      dpi = 600
    )
  }
)
