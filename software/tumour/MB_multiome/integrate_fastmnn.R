# ==============================================================================
# Integrate Taylor lab MB tumour scRNA-seq samples using fastMNN.
# ==============================================================================

library(argparse)
library(tidyverse)
library(patchwork)
library(Seurat)
library(SeuratWrappers)
library(Signac)

# parse command line arguments
parser <- ArgumentParser()

parser$add_argument(
  # input directory with the mtx/barcodes/genes
  "--in_dir",
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
    "--in_dir", "/.mounts/labs/pailab/private/projects/MB_multiome",
    "--out_dir", "/.mounts/labs/pailab/private/icheong/CBL_scRNAseq/results/tumour/MB_multiome/20250206"
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

#' Integrate a Seurat object with/without fastMNN. Copied and modified from
#' `test_integ_methods.R`.
#' 
#' @param srat_list List of Seurat objects to integrate.
#' @param fastMNN If `TRUE` (default), integrate Seurat objects with fastMNN.
#'   Otherwise, combine the Seurat objects without integrating (using `merge`).
#' @param features Number of variable features to use for integration or list of
#'   features.
#' @param npcs Number of principal components to calculate in `RunPCA`.
#' @param ndims Number of principal components to use for downstream analysis.
#'
#' @return Seurat object after fastMNN integration.
#' 
integ_srat <- function(
  srat_list,
  fastMNN = TRUE,
  features = 3000,
  npcs = 100,
  ndims = 50
) {
  if (fastMNN) {
    # run fastMNN
    message(sprintf(
      "***Integrating %s Seurat objects with `fastMNN`***",
      length(srat_list)
    ))
    integ_srat <- SeuratWrappers::RunFastMNN(
      srat_list,
      assay = "SCT",
      features = features
    )
    reduction <- "mnn"
  } else {
    # merge instead of integrate
    message(sprintf(
      "***Merging %s Seurat objects with `merge`; no integration being performed***",
      length(srat_list)
    ))
    integ_srat <- merge(x = srat_list[[1]], y = srat_list[-1])

    # set variable features manually
    if (is.numeric(features)) {
      features <- SelectIntegrationFeatures(srat_list, nfeatures = features)
    }
    VariableFeatures(integ_srat) <- features

    # re-calculate residuals for scale.data slot
    integ_srat <- GetResidual(integ_srat, features = features)

    # run PCA (required for UMAP)
    integ_srat <- RunPCA(integ_srat, npcs = npcs)
    reduction <- "pca"
  }

  # UMAP
  integ_srat <- RunUMAP(integ_srat, reduction = reduction, dims = 1:ndims)

  # clustering
  integ_srat <- FindNeighbors(integ_srat, reduction = reduction, dims = 1:ndims) %>% 
    FindClusters(.)

  return(integ_srat)
}

source("/.mounts/labs/pailab/private/icheong/CBL_scRNAseq/software/utilities/plotting.R")

make_plots <- function(
  srat,
  path
) {
  # original samples
  umap_sample <- DimPlot(
    srat,
    reduction = "umap",
    group.by = "sample_id",
    label = TRUE,
    label.size = 2,
    repel = TRUE,
    raster = FALSE
  ) + 
    NoLegend()
  ggsave(
    plot = umap_sample,
    filename = "umap_sample.png",
    path = path,
    width = 8,
    height = 8,
    units = "in",
    dpi = 600
  )

  # highlight each of the samples on the plot
  umap_sample_hl <- highlight_DimPlot(
    srat,
    highlight_by = "sample_id",
    reduction = "umap",
    sizes.highlight = 0.01,
    label = FALSE,
    raster = FALSE
  )
  nc <- ceiling(sqrt(length(umap_sample_hl)))
  nr <- ceiling(length(umap_sample_hl) / nc)
  ggsave(
    plot = wrap_plots(umap_sample_hl, ncol = nc),
    filename = "umap_sample_split.png",
    path = path,
    width = nc * 4,
    height = nr * 4,
    units = "in",
    dpi = 300
  )

  # subtypes
  umap_subtype <- DimPlot(
    srat,
    reduction = "umap",
    group.by = "subtype",
    label = TRUE,
    repel = TRUE,
    raster = FALSE
  ) + 
    NoLegend()
  ggsave(
    plot = umap_subtype,
    filename = "umap_subtype.png",
    path = path,
    width = 8,
    height = 8,
    units = "in",
    dpi = 600
  )

  # highlight each of the subtypes on the plot
  umap_subtype_hl <- highlight_DimPlot(
    srat,
    highlight_by = "subtype",
    reduction = "umap",
    sizes.highlight = 0.01,
    label = FALSE,
    raster = FALSE
  )
  nc <- ceiling(sqrt(length(umap_subtype_hl)))
  nr <- ceiling(length(umap_subtype_hl) / nc)
  ggsave(
    plot = wrap_plots(umap_subtype_hl, ncol = nc),
    filename = "umap_subtype_split.png",
    path = path,
    width = nc * 4,
    height = nr * 4,
    units = "in",
    dpi = 300
  )

  # clusters
  umap_clusters <- DimPlot(
    srat,
    reduction = "umap",
    group.by = "seurat_clusters",
    label = TRUE,
    repel = TRUE,
    raster = FALSE
  ) + 
    NoLegend()
  ggsave(
    plot = umap_clusters,
    filename = "umap_clusters.png",
    path = path,
    width = 8,
    height = 8,
    units = "in",
    dpi = 600
  )

  # highlight each of the clusters on the plot
  umap_clusters_hl <- highlight_DimPlot(
    srat,
    highlight_by = "seurat_clusters",
    reduction = "umap",
    sizes.highlight = 0.01,
    label = FALSE,
    raster = FALSE
  )
  nc <- ceiling(sqrt(length(umap_clusters_hl)))
  nr <- ceiling(length(umap_clusters_hl) / nc)
  ggsave(
    plot = wrap_plots(umap_clusters_hl, ncol = nc),
    filename = "umap_clusters_split.png",
    path = path,
    width = nc * 4,
    height = nr * 4,
    units = "in",
    dpi = 300
  )

  # combine umaps on one plot
  .plt <- wrap_plots(umap_sample, umap_subtype, umap_clusters, nrow = 1) & 
    NoAxes()
  ggsave(
    plot = .plt,
    filename = "umaps.png",
    path = path,
    width = 24,
    height = 8,
    units = "in",
    dpi = 600
  )

  # bar plots (from `cluster_barplot.R`)
  .plt1 <- cluster_barplot(
    srat,
    split.by = "seurat_clusters",
    group.by = "sample_id",
    position = "fill",
    width = 1
  ) + 
    theme(axis.text.x = element_blank())
  .plt2 <- cluster_barplot(
    srat,
    split.by = "subtype",
    group.by = "seurat_clusters"
  )
  .plt <- wrap_plots(.plt1, .plt2, nrow = 2)
  ggsave(
    plot = .plt,
    filename = "stacked_barplots.png",
    path = path,
    width = 15,
    height = 10,
    units = "in",
    dpi = 600
  )
}

# ------------------------------------------------------------------------------
# load datasets

# JSON file from Gabrielle Persad (Stein lab) with list of samples to use
qc_filters <- jsonlite::fromJSON(file.path(arg_list$in_dir, "MB_qc_filters.json")) %>%
  pluck(1)

# load tumour samples and prep for integration
srat_list <- map(
  .x = names(qc_filters),
  .f = \(sample_id) {
    message(sprintf("***Reading in %s***", sample_id))
    srat <- readRDS(file.path(arg_list$in_dir, "input", paste0(sample_id, ".rds")))

    # set SCT as default (active) assay
    DefaultAssay(srat) <- "SCT"

    # reduce size of object
    srat <- DietSeurat(
      srat,
      assays = c("RNA", "SCT")
    )

    # add sample ID to cell names to ensure uniqueness when integrating
    srat <- RenameCells(
      srat,
      add.cell.id = sample_id
    )

    # add subtype to metadata
    srat$subtype <- qc_filters[[sample_id]]$subtype

    return(srat)
  }
)

# ------------------------------------------------------------------------------
# integrate datasets with/without fastMNN

# integrate with fastMNN and save output
mb_mnn <- integ_srat(
  srat_list = srat_list,
  fastMNN = TRUE,
  features = 3000,
  npcs = 100,
  ndims = 25
)

message("***Saving integrated object as mb_fastmnn.qs***")
qs::qsave(
  x = mb_mnn,
  file = file.path(arg_list$out_dir, "fastmnn", "mb_fastmnn.qs")
)

# merge datasets only (no fastMNN integration) and save output
mb_mrg <- integ_srat(
  srat_list = srat_list,
  fastMNN = FALSE,
  features = 3000,
  npcs = 100,
  ndims = 25
)

message("***Saving merged object as mb_merged.qs***")
qs::qsave(
  x = mb_mrg,
  file = file.path(arg_list$out_dir, "merged", "mb_merged.qs")
)

# ------------------------------------------------------------------------------
# plot datasets

# plots for fastMNN samples
make_plots(srat = mb_mnn, path = file.path(arg_list$out_dir, "fastmnn"))

# plots for merged samples
make_plots(srat = mb_mrg, path = file.path(arg_list$out_dir, "merged"))

# ------------------------------------------------------------------------------

message("\n***SESSION INFO***\n")
print(sessionInfo())
