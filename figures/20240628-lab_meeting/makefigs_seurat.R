# ==============================================================================
# These figures were generated for the lab meeting on 2024-06-28. This script
# should be run in the `scrnaseq_env` conda environment.
# ==============================================================================

library(tidyverse)
library(patchwork)
library(ggrepel)
library(Seurat)

# >>> load Seurat >>>

# Seurat objects for import
srat_names <- c("full_cca", "rl_cca", "ubc_cca")
srat_qs <- c(
  # full cerebellum CCA integration
  "/.mounts/labs/pailab/private/llau/results/integrated/20240516/cca/20240516_cca_integ.qs",
  # RL lineage CCA integration
  "/.mounts/labs/pailab/private/llau/results/integrated/20240524/25_pc_without_luo/25_pc_rl.qs",
  # UBC CCA integration
  "/.mounts/labs/pailab/private/llau/results/integrated/20240527/ubc_subset.qs"
)

# import Seurat objects
all_srat <- map(
  .x = srat_qs,
  .f = qs::qread
)
names(all_srat) <- srat_names
print(all_srat)

# set species levels and default assay to "RNA"
all_srat <- map(
  .x = all_srat,
  .f = \(x) {
    x$species <- factor(x$species, levels = c("mouse", "human"))
    DefaultAssay(x) <- "RNA"
    return(x)
  }
)

# add UBC subclusters to RL integration metadata
md <- merge(
  x = all_srat$rl_cca[[]],
  y = all_srat$ubc_cca[[]] %>%
    select(snn_res.0.3) %>%
    rename(ubc_subclusters = snn_res.0.3),
  by = "row.names",
  all.x = TRUE
) %>%
  column_to_rownames(var = "Row.names")
all_srat$rl_cca@meta.data <- md[rownames(all_srat$rl_cca[[]]), ]

# <<< load Seurat <<<

# colour palette for species
species_cols <- c("grey", scales::pal_hue()(1))

# import functions
source("/.mounts/labs/pailab/private/icheong/CBL_scRNAseq/software/utilities/plotting.R")

# ------------------------------------------------------------------------------
# UMAPs of RL lineage integration

rl_cluster <- "snn_res.0.4"
rl_cluster_cols <- pals::cols25(length(table(all_srat$rl_cca[[rl_cluster]])))

# coloured by cluster
.plt <- DimPlot(
  all_srat$rl_cca,
  reduction = "umap",
  group.by = rl_cluster,
  cols = rl_cluster_cols,
  label = TRUE,
  repel = TRUE
) + 
  NoLegend() + 
  labs(title = NULL)
ggsave(
  file = "rl_cca_cluster.png",
  plot = .plt,
  width = 5.5,
  height = 5,
  units = "in",
  dpi = 600
)

# coloured by species
.plt <- DimPlot(
  all_srat$rl_cca,
  reduction = "umap",
  cells.highlight = WhichCells(all_srat$rl_cca, expression = species == "human"),
  sizes.highlight = 0.01,
  label = FALSE
) + 
  labs(title = NULL) + 
  scale_colour_manual(
    labels = c("mouse", "human"),
    values = species_cols
  )
ggsave(
  file = "rl_cca_species.png",
  plot = .plt,
  width = 6,
  height = 5,
  units = "in",
  dpi = 600
)

# coloured by UBC subclusters
.plt <- DimPlot(
  all_srat$rl_cca,
  reduction = "umap",
  group.by = "ubc_subclusters",
  label = TRUE,
  repel = TRUE,
  cols = scales::pal_hue()(length(table(all_srat$rl_cca$ubc_subclusters))),
  na.value = "grey"
) + 
  NoLegend() + 
  labs(title = NULL)
ggsave(
  file = "rl_ubc_subclusters.png",
  plot = .plt,
  width = 4.5,
  height = 4,
  units = "in",
  dpi = 600
)

# plot EOMES and LMX1A expression
walk(
  .x = c("EOMES", "LMX1A"),
  .f = \(gene) {
    .plt <- FeaturePlot(
      all_srat$rl_cca,
      features = gene,
      order = TRUE,
      min.cutoff = "q10",
      max.cutoff = "q90"
    )
    ggsave(
      filename = paste0(gene, ".png"),
      plot = .plt,
      width = 5,
      height = 4,
      units = "in",
      dpi = 600
    )
  }
)

# ------------------------------------------------------------------------------
# bar plots

# proportion of cells in each RL cluster
.plt <- cluster_barplot(
  object = all_srat$rl_cca,
  split.by = rl_cluster,
  group.by = "species",
  position = "fill"
) + 
  scale_fill_manual(values = rl_cluster_cols) + 
  theme(legend.title = element_blank())
ggsave(
  filename = "rl_cluster_bar.png",
  plot = .plt,
  width = 4,
  height = 5,
  units = "in",
  dpi = 600
)

# proportion of cells in each UBC subcluster
.plt <- cluster_barplot(
  object = all_srat$ubc_cca,
  split.by = "snn_res.0.3",
  group.by = "species",
  position = "fill"
) + 
  theme(legend.title = element_blank())
ggsave(
  filename = "ubc_cluster_bar.png",
  plot = .plt,
  width = 3,
  height = 5,
  units = "in",
  dpi = 600
)
