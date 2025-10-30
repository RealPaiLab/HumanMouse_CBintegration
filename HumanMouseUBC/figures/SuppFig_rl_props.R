# ==============================================================================
# Supplementary figure of manuscript showing the proportion tests for the RL
# lineage cell types. This script should be run using the scrnaseq_env conda
# environment.
# ==============================================================================

library(tidyverse)
library(patchwork)
library(Seurat)

source("./utils.R")
source("../../software/utilities/cell_labelling.R")
source("../../software/utilities/plotting.R")
source("../../software/utilities/propeller_helpers.R")

out_dir <- file.path(
  "/.mounts/labs/pailab/private/projects/HumanMouseUBC/figures/SuppFig_rl_props",
  format(Sys.Date(), "%Y%m%d")
)
if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

# colour palette
my_pals <- get_custom_pals()

# ------------------------------------------------------------------------------
# proportion of UBCs in RL

# load Seurat object for RL lineage cells and isolated UBCs
srat_qs <- get_srat_paths()
rl_srat <- load_srat(srat_qs["rl"]) %>%
  pluck("rl")

# get cell type annotations (from `cell_labelling.R`)
rl_srat <- label_rl_lineage_integration(rl_srat)

# collapse UBC subclusters into one cluster, set species factor levels
rl_srat@meta.data <- mutate(
  rl_srat[[]],
  broad_annot = case_when(
    str_starts(cell_type_annot, "UBC") ~ "UBC",
    .default = cell_type_annot
  ),
  species = factor(species, levels = c("mouse", "human"))
)

# filter out control cells (endothelial cells, microglia, oligodendrocytes)
rl_srat <- subset(
  x = rl_srat,
  subset = broad_annot %in% c("endothelial", "microglia", "oligodendrocyte/OPC"),
  invert = TRUE
)

# proportion of UBCs in human vs mouse RL
prop_bar <- cluster_barplot(
  object = rl_srat,
  split.by = "broad_annot",
  group.by = "species",
  position = "fill"
) + 
  labs(fill = "general\ncell type") + 
  scale_fill_manual(values = my_pals$rl_integ_annot) + 
  theme_classic2()
ggsave(
  filename = "cluster_prop_bar.pdf",
  plot = prop_bar,
  path = out_dir,
  width = 4,
  height = 4,
  units = "in"
)

# ------------------------------------------------------------------------------
# scProportionTest

# load scProportionTest results
scpt_res <- read_csv("/.mounts/labs/pailab/private/icheong/CBL_scRNAseq/results/integrated/20240822/rl_broad_annot_results.csv")

# add significance column and set plot order
log2FD_threshold <- log2(2)
FDR_threshold <- 0.05
scpt_res <- mutate(
  scpt_res,
  significance = case_when(
    FDR < FDR_threshold & obs_log2FD > log2FD_threshold ~ "human",
    FDR < FDR_threshold & obs_log2FD < -log2FD_threshold ~ "mouse",
    .default = "n.s."
  ),
  significance = factor(significance, levels = c("mouse", "n.s.", "human")),
  clusters = fct_reorder(clusters, obs_log2FD, .desc = FALSE)
)

# plot scProportionTest results
scpt_plt <- ggplot(
  data = scpt_res,
  mapping = aes(
    x = obs_log2FD,
    y = clusters,
    xmin = boot_CI_2.5,
    xmax = boot_CI_97.5,
    colour = significance
  )
) + 
  geom_pointrange(size = 0.25) + 
  geom_vline(
    xintercept = c(log2FD_threshold, -log2FD_threshold),
    lty = "dashed"
  ) + 
  geom_vline(xintercept = 0, lty = "solid") + 
  labs(x = "observed log2(FD)", y = "cell type", colour = "enriched") + 
  scale_colour_manual(values = c(my_pals$species[1], "black", my_pals$species[2])) + 
  theme_classic2() + 
  theme(legend.position = "bottom")
ggsave(
  filename = "permutation_plot.pdf",
  plot = scpt_plt,
  path = out_dir,
  width = 5,
  height = 3,
  units = "in"
)

# ------------------------------------------------------------------------------
# propeller

# set same cluster plotting order as scProportionTest
cluster_levels <- levels(scpt_res$clusters)

# load `propeller` sample proportions
propeller_props <- read_csv("/.mounts/labs/pailab/private/icheong/CBL_scRNAseq/results/integrated/20240823/rl_proportions.csv") %>%
  mutate(
    dataset_name = str_replace(
      string = dataset_name,
      pattern = "_([:alpha:]+)",
      replacement = " (\\1)"
    ),
    species = factor(species, levels = c("mouse", "human")),
    clusters = factor(clusters, levels = cluster_levels)
  )

# load propeller
propeller_res <- read_csv("/.mounts/labs/pailab/private/icheong/CBL_scRNAseq/results/integrated/20240823/rl_results.csv")
fdrs <- get_fdrs_for_plotting(propeller_res) %>%
  mutate(clusters = factor(clusters, levels = cluster_levels))

# plot `propeller` results for UBCs
propeller_plt <- boxplot_props_per_cluster(
  props = propeller_props,
  facet_var = clusters,
  facet_ncol = 4,
  fdrs = fdrs,
  fdr_label_size = 8.8
) + 
  labs(y = "proportion of all\nRL lineage cells", shape = "dataset") + 
  theme_classic2()
ggsave(
  filename = "propeller_plot.pdf",
  plot = propeller_plt,
  path = out_dir,
  width = 10,
  height = 4,
  units = "in"
)

# ------------------------------------------------------------------------------
# final figure

layout <- c(
  area(t = 1, l = 1, r = 2),
  area(t = 1, l = 3, r = 5),
  area(t = 2, l = 1, r = 5)
)

.plt <- wrap_plots(
  free(prop_bar),
  free(scpt_plt),
  free(propeller_plt)
) + 
  plot_layout(design = layout)

walk(
  .x = c("png", "pdf"),
  .f = \(device) {
    ggsave(
      filename = paste0("combined.", device),
      plot = .plt,
      path = out_dir,
      width = 9,
      height = 6,
      units = "in",
      dpi = 600
    )
  }
)
