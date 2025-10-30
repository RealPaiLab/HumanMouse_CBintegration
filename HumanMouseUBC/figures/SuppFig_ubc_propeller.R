# ==============================================================================
# Supplementary figure of manuscript showing the propeller proportion test
# results. This script should be run using the scrnaseq_env conda environment.
# ==============================================================================

library(tidyverse)
library(Seurat)

source("./utils.R")
source("../../software/utilities/propeller_helpers.R")

out_dir <- file.path(
  "/.mounts/labs/pailab/private/projects/HumanMouseUBC/figures/SuppFig_ubc_propeller",
  format(Sys.Date(), "%Y%m%d")
)
if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

# ------------------------------------------------------------------------------
# propeller

# load `propeller` sample proportions
propeller_props <- read_csv("/.mounts/labs/pailab/private/icheong/CBL_scRNAseq/results/integrated/20240606/four_dataset_cca_ubc_proportions.csv") %>%
  mutate(
    dataset_name = str_replace(
      string = dataset_name,
      pattern = "_([:alpha:]+)",
      replacement = " (\\1)"
    ),
    species = factor(species, levels = c("mouse", "human")),
    clusters = paste0("iUBC", clusters)
  )

# load propeller
propeller_res <- read_csv("/.mounts/labs/pailab/private/icheong/CBL_scRNAseq/results/integrated/20240606/four_dataset_cca_ubc_results.csv")
fdrs <- get_fdrs_for_plotting(propeller_res) %>%
  mutate(
    clusters = paste0("iUBC", clusters)
  )

# plot `propeller` results for UBCs
propeller_plt <- boxplot_props_per_cluster(
  props = propeller_props,
  facet_var = clusters,
  facet_ncol = 3,
  fdrs = fdrs,
  fdr_label_size = 8.8
) + 
  labs(y = "proportion of all RL lineage cells", shape = "dataset") + 
  theme_classic2()

walk(
  .x = c("png", "pdf"),
  .f = \(device) {
    ggsave(
      filename = paste0("ubc_propeller_plot.", device),
      plot = propeller_plt,
      path = out_dir,
      width = 6,
      height = 5,
      units = "in",
      dpi = 600
    )
  }
)
