# ==============================================================================
# Supplementary figure of manuscript showing UMAPs of the RL lineage cells and
# RL cell type markers (e.g. DotPlot). This script should be run using the
# scrnaseq_env conda environment.
# ==============================================================================

library(tidyverse)
library(patchwork)
library(Seurat)

source("./utils.R")
source("../../software/utilities/cell_labelling.R")
source("../../software/utilities/plotting.R")

out_dir <- file.path(
  "/.mounts/labs/pailab/private/projects/HumanMouseUBC/figures/SuppFig_rl_markers",
  format(Sys.Date(), "%Y%m%d")
)
if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

# colour palette
my_pals <- get_custom_pals()

# ------------------------------------------------------------------------------
# dot plot of cell type markers

# load Seurat object for integrated dataset
srat_qs <- get_srat_paths()
integ_srat <- load_srat(srat_qs["rl"]) %>%
  pluck("rl")

# get cell type annotations (from `cell_labelling.R`)
integ_srat <- label_rl_lineage_integration(integ_srat)

# collapse UBC subclusters into one cluster
integ_srat@meta.data <- mutate(
  integ_srat[[]],
  broad_annot = case_when(
    str_starts(cell_type_annot, "UBC") ~ "UBC",
    .default = cell_type_annot
  ),
  broad_annot = factor(
    broad_annot,
    levels = c("RL", "UBC", "GCP", "GC", "oligodendrocyte/OPC", "microglia", "endothelial")
  ),
  dataset_name = str_remove(
    string = dataset_name,
    pattern = "full_cerebellum_"
  )
)

# normalize the RNA assay (the assay looks like it should be normalized, but
# we'll just re-run it anyways to be sure)
DefaultAssay(integ_srat) <- "RNA"
integ_srat <- NormalizeData(integ_srat)

# marker genes for the general cell types; selected from /.mounts/labs/pailab/private/icheong/CBL_scRNAseq/scrnaseq_Leo/utilities/cell_gene_mapping.csv
cell_type_markers <- c(
  "MKI67", # RL, GCP
  "WLS", # RL, UBC progenitors
  "EOMES", "LMX1A", "OTX2", # UBC
  "RBFOX3", "ATOH1", # GCP
  "PAX6", # differentiated cells (UBC/GC)
  "RELN", # GC
  "PDGFRA", "TNR", # OPC/oligodendrocytes
  "TREM2", # microglia
  "CLDN5", "ITM2A" # endothelial
)

# dot plot
dot_plot <- DotPlot(
  integ_srat,
  assay = "RNA",
  features = rev(cell_type_markers),
  group.by = "broad_annot"
) + 
  scale_colour_gradient(
    low = "lightgrey",
    high = "blue",
    breaks = \(x) {quantile(x, c(0.1, 0.9))}, # use 10/90 percentile as breaks
    labels = c("low", "high")
  ) + 
  coord_flip() + 
  theme_classic2() + 
  theme(
    axis.text.x = element_text(hjust = 1, angle = 30),
    legend.ticks = element_blank()
  )
ggsave(
  filename = "dot_plot.pdf",
  plot = dot_plot,
  path = out_dir,
  width = 5,
  height = 6,
  units = "in"
)
