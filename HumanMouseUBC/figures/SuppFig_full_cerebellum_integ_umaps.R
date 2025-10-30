# ==============================================================================
# Supplementary figure of manuscript. This script should be run using the
# scrnaseq_env conda environment.
# ==============================================================================

library(tidyverse)
library(patchwork)
library(cowplot)
library(Seurat)

source("./utils.R")

out_dir <- file.path(
  "/.mounts/labs/pailab/private/projects/HumanMouseUBC/figures/SuppFig_full_cerebellum_integ_umaps",
  format(Sys.Date(), "%Y%m%d")
)
my_pals <- get_custom_pals()

# ------------------------------------------------------------------------------
# UMAPs of individual datasets

# load Seurat objects for individual datasets
srat_qs <- get_srat_paths()
srat_list <- load_srat(srat_qs[c("aldinger_human", "sepp_human", "vladoiu_mouse", "sepp_mouse")])

# get all common cell types and set as factor levels
common_cell_types <- map(
  .x = srat_list,
  .f = \(srat) {unique(srat$common_cell_name)}
) %>%
  list_c() %>%
  unique() %>%
  purrr::discard(.p = is.na) %>%
  sort()

srat_list <- map(
  .x = srat_list,
  .f = \(srat) {
    srat$common_cell_name <- factor(
      srat$common_cell_name,
      levels = common_cell_types
    )
    return(srat)
  }
)

# plot UMAPs for each individual dataset
indiv_umaps <- map2(
  .x = srat_list,
  .y = names(srat_list),
  .f = \(srat, dataset) {
    plt_title <- str_to_sentence(dataset) %>%
      str_replace(pattern = "_([:alpha:]+)", replacement = " (\\1)")

    .plt <- DimPlot(
      srat,
      reduction = "umap",
      group.by = "common_cell_name",
      label = TRUE,
      label.size = 3,
      repel = TRUE,
      raster = FALSE
    ) + 
      labs(title = plt_title, x = "UMAP 1", y = "UMAP 2") + 
      scale_colour_manual(
        drop = FALSE,
        values = scales::pal_hue()(25)
      )
    
    # hacky method so legends can be combined in downstream patchwork (normally
    # you'd set `geom_point(show.legend = TRUE)`, but this not available with
    # `DimPlot`); also need to deal with NA
    .plt[[1]]$layers[[1]]$show.legend <- TRUE
    .plt[[1]]$data$common_cell_name <- addNA(.plt[[1]]$data$common_cell_name)

    ggsave(
      filename = paste0(dataset, "_umap.pdf"),
      plot = .plt,
      path = out_dir,
      width = 10,
      height = 5,
      units = "in"
    )

    return(.plt)
  }
)

# free up some RAM
rm(srat_list)

# ------------------------------------------------------------------------------
# UMAP of full cerebellum integration

# load Seurat object for integrated dataset
integ_srat <- load_srat(srat_qs["full"]) %>%
  pluck("full")
integ_srat$dataset_name <- str_remove(
  string = integ_srat$dataset_name,
  pattern = "full_cerebellum_"
)

# plot UMAP by cell type
full_umap <- DimPlot(
  integ_srat,
  reduction = "umap",
  group.by = "common_cell_name",
  label = TRUE,
  label.size = 3,
  repel = TRUE,
  raster = FALSE
) + 
  labs(title = "Human+mouse integration", x = "UMAP 1", y = "UMAP 2") + 
  guides(colour = guide_legend(
    theme = theme(
      legend.text = element_text(size = rel(0.6), margin = margin(l = unit(1, "points"))),
      # legend.key.spacing.x = unit(6, "points"),
      legend.key.spacing.y = unit(2, "points")
    ),
    position = "bottom",
    override.aes = list(size = 3),
    ncol = 2
  )) + 
  theme_classic2()
ggsave(
  filename = "integ_full_cell_type.pdf",
  plot = full_umap,
  path = out_dir,
  width = 5,
  height = 7.5,
  units = "in"
)

# plot UMAP by species
full_umap_species <- DimPlot(
  integ_srat,
  reduction = "umap",
  cells.highlight = WhichCells(integ_srat, expression = species == "human"),
  sizes.highlight = 0.01,
  label = FALSE,
  raster = FALSE
) + 
  labs(title = "Human+mouse (species)", x = "UMAP 1", y = "UMAP 2") + 
  scale_colour_manual(
    labels = c("mouse", "human"),
    values = my_pals$species
  ) + 
  theme_classic2() + 
  theme(
    legend.position = "bottom",
    legend.justification = c(0.5, 1),
    legend.key.spacing.x = unit(8, "points")
  )
ggsave(
  filename = "integ_full_species.pdf",
  plot = full_umap_species,
  path = out_dir,
  width = 5,
  height = 5.5,
  units = "in"
)

# plot UMAP by dataset
full_umap_dataset <- DimPlot(
  integ_srat,
  reduction = "umap",
  group.by = "dataset_name",
  label = FALSE,
  raster = FALSE
) + 
  labs(title = "Human+mouse (dataset)", x = "UMAP 1", y = "UMAP 2") + 
  scale_colour_manual(
    values = alpha(scales::pal_hue()(4), alpha = 0.2)
  ) + 
  guides(colour = guide_legend(
    override.aes = list(size = 3, alpha = 1),
    nrow = 2
  )) + 
  theme_classic2() + 
  theme(
    legend.position = "bottom",
    legend.justification = c(0.5, 1),
    legend.key.spacing.x = unit(12, "points")
  )
ggsave(
  filename = "integ_full_dataset.pdf",
  plot = full_umap_dataset,
  path = out_dir,
  width = 5,
  height = 5.5,
  units = "in"
)

# ------------------------------------------------------------------------------
# final figure

indiv_umaps <- map(
  .x = indiv_umaps,
  .f = \(plt) {
    plt <- plt + 
      theme_classic2() + 
      theme(legend.position = "none")
  }
)

.plt <- plot_grid(
  plot_grid(
    NULL,
    indiv_umaps$aldinger_human,
    indiv_umaps$sepp_human,
    NULL,
    rel_widths = c(1, 2, 2, 1),
    ncol = 4
  ),
  plot_grid(
    NULL,
    indiv_umaps$vladoiu_mouse,
    indiv_umaps$sepp_mouse,
    NULL,
    rel_widths = c(1, 2, 2, 1),
    ncol = 4
  ),
  plot_grid(
    full_umap,
    full_umap_species,
    full_umap_dataset,
    align = "h",
    axis = "tb",
    ncol = 3
  ),
  nrow = 3,
  rel_heights = c(0.26, 0.26, 0.48)
)

walk(
  .x = c("png", "pdf"),
  .f = \(device) {
    ggsave(
      filename = paste("combined", device, sep = "."),
      plot = .plt,
      path = out_dir,
      width = 10.5,
      height = 12.25,
      units = "in",
      dpi = 600
    )
  }
)
