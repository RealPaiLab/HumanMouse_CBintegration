# ==============================================================================
# Supplementary figure of manuscript showing the `clustree` visualization from
# clustering the UBCs at different resolutions. This script should be run using
# the clustree conda environment.
# ==============================================================================

library(tidyverse)
library(Seurat)
library(clustree)

source("./utils.R")

out_dir <- file.path(
  "/.mounts/labs/pailab/private/projects/HumanMouseUBC/figures/SuppFig_ubc_clustree",
  format(Sys.Date(), "%Y%m%d")
)
if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

# ------------------------------------------------------------------------------
# clustering tree plot

# load Seurat object for integrated dataset
srat_qs <- get_srat_paths()
ubc_srat <- load_srat(srat_qs["ubc"]) %>%
  pluck("ubc")

.plt <- clustree(ubc_srat, prefix = "snn_res.") + 
  labs(title = "Clustering tree of UBCs") + 
  guides(
    colour = guide_legend(title = "cluster\nresolution", order = 1),
    edge_alpha = guide_legend(title = "in-proportion", order = 2),
    edge_colour = guide_legend(title = "number\nof cells", order = 3),
    size = guide_legend(title = "number\nof cells", order = 4)
  ) + 
  theme(
    plot.title = element_text(face = "plain", hjust = 0.5),
    plot.title.position = "plot",
    plot.margin = margin(5.5, 5.5, 5.5, 5.5, "points")
  )

walk(
  .x = c("png", "pdf"),
  .f = \(device) {
    ggsave(
      filename = paste0("clustree.", device),
      plot = .plt,
      path = out_dir,
      width = 8,
      height = 9,
      units = "in",
      dpi = 600
    )
  }
)
