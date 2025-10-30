# ==============================================================================
# Test different clustering resolutions to determine the number of tumour
# clusters.
# ==============================================================================

library(argparse)
library(tidyverse)
library(patchwork)
library(Seurat)
library(clustree) # <-- USE CLUSTREE CONDA ENVIRONMENT

# parse command line arguments
parser <- ArgumentParser()

parser$add_argument(
  # input directory with the mtx/barcodes/genes
  "--srat_file",
  default = NULL,
  required = TRUE
)
parser$add_argument(
  # output directory
  "--out_dir",
  default = NULL,
  required = TRUE
)

if (interactive()) {
  # for testing and troubleshooting
  arg_list <- parser$parse_args(c(
    "--srat_file", "/.mounts/labs/pailab/private/icheong/CBL_scRNAseq/results/tumour/Vladoiu/20230510/mb_mnn.rds",
    "--out_dir", "/.mounts/labs/pailab/private/icheong/CBL_scRNAseq/results/tumour/Vladoiu/20250211"
  ))
} else {
  arg_list <- parser$parse_args()
}

message(sprintf("Saving files to %s", getwd()))
if (!dir.exists(arg_list$out_dir)) {
  dir.create(arg_list$out_dir, recursive = TRUE)
}

# ------------------------------------------------------------------------------
# functions

source("/.mounts/labs/pailab/private/icheong/CBL_scRNAseq/software/utilities/plotting.R")

# ------------------------------------------------------------------------------
# load dataset

srat <- readRDS(arg_list$srat_file)

# remove previously determined clusters
srat@meta.data <- srat[[]] %>%
  select(!contains("snn_res"))

# turn subtype into factor
srat$subtype <- factor(srat$subtype, levels = c("SHH", "G3", "G4"))

# ------------------------------------------------------------------------------
# cluster at different resolutions

# resolutions for clustering
test_res <- c(
  seq(from = 0.1, to = 0.5, by = 0.05),
  seq(from = 0.6, to = 1.2, by = 0.1)
)
for (res in test_res) {
  message(sprintf("***Finding clusters at resolution: %s***", res))
  srat <- FindClusters(srat, resolution = res)
}

# save Seurat object and the metadata
qs::qsave(x = srat, file = file.path(arg_list$out_dir, "mb_mnn.qs"))
write.table(
  x = srat[[]],
  file = file.path(arg_list$out_dir, "mb_mnn_metadata.tsv"),
  quote = FALSE,
  sep = "\t"
)

# plot clustree
.plt <- clustree(
  x = srat,
  prefix = "SCT_snn_res.",
  prop_filter = 0.05
)
ggsave(
  plot = .plt,
  filename = "clustree.png",
  path = arg_list$out_dir,
  width = 10,
  height = 15,
  units = "in",
  dpi = 600
)

# plot UMAPs and bar plots for each resolution
walk(
  .x = test_res,
  .f = \(res) {
    highlight_by <- sprintf("SCT_snn_res.%s", res)

    # UMAPs highlighting each cluster
    .plt <- highlight_DimPlot(
      srat,
      highlight_by = highlight_by,
      reduction = "umap",
      sizes.highlight = 0.01,
      label = FALSE,
      raster = FALSE
    )
    nc <- ceiling(sqrt(length(.plt)))
    nr <- ceiling(length(.plt) / nc)
    ggsave(
      plot = wrap_plots(.plt, ncol = nc),
      filename = sprintf("umap_%.02f_split.png", res),
      path = arg_list$out_dir,
      width = nc * 4,
      height = nr * 4,
      units = "in",
      dpi = 300
    )

    # UMAP with all clusters
    .plt <- DimPlot(
      srat,
      reduction = "umap",
      group.by = highlight_by,
      label = TRUE,
      repel = TRUE
    ) + 
      NoLegend()
    ggsave(
      plot = .plt,
      filename = sprintf("umap_%.02f.png", res),
      path = arg_list$out_dir,
      width = 5,
      height = 5,
      units = "in",
      dpi = 600
    )

    # bar plot to show subtype in each cluster
    .plt <- cluster_barplot(
      srat,
      split.by = "subtype",
      group.by = highlight_by,
      position = "stack",
      width = 0.75
    )
    ggsave(
      plot = .plt,
      filename = sprintf("bar_subtype_%.02f.png", res),
      path = arg_list$out_dir,
      width = (length(table(srat[[highlight_by]])) * 0.4) + 1,
      height = 5,
      units = "in",
      dpi = 600
    )
  }
)

# ------------------------------------------------------------------------------

message("\n***SESSION INFO***\n")
print(sessionInfo())
