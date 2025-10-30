# ==============================================================================
# Test different integration methods on the data.
# ==============================================================================

# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(argparse)
library(tidyverse)
library(patchwork)
library(Seurat)
library(harmony)

# parse command line arguments
parser <- ArgumentParser()

parser$add_argument(
  # input directory with the mtx/barcodes/genes
  "--srat_rds",
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
  args <- parser$parse_args(c(
    "--srat_rds", "/CBL_scRNAseq/results/tumour/Vladoiu/20230508/vladoiu_mb_merge.rds",
    "--out_dir", "/CBL_scRNAseq/results/tumour/Vladoiu/20230510/"
  ))
} else {
  args <- parser$parse_args()
}

message(sprintf("Saving files to %s", getwd()))
if (!dir.exists(args$out_dir)) {
  dir.create(args$out_dir, recursive = TRUE)
}

# ------------------------------------------------------------------------------
# load Seurat

srat <- readRDS(args$srat_rds)

srat$subtype <- stringr::str_replace(
  string = srat$orig.ident,
  pattern = "^[:alnum:]*_",
  replacement = ""
)

# ------------------------------------------------------------------------------
# functions


#' Use Seurat CCA to integrate a Seurat object.
#'
#' @param srat Seurat object.
#' @param out_dir Directory to save output plots.
#' @param integ_var Metadata column to integrate on.
#' @param npcs Number of principal components to calculate in `RunPCA`.
#' @param ndims Number of principal components to use for downstream analysis.
#'
#' @return Seurat object after CCA integration.
#'
integ_cca <- function(
  srat,
  out_dir,
  integ_var = "orig.ident",
  npcs = 100,
  ndims = 50
) {
  # run sctransform
  srat <- SCTransform(
    srat,
    assay = "RNA",
    variable.features.n = 5000,
    vars.to.regress = "CC.Difference",
    return.only.var.genes = FALSE
  )
  
  srat_list <- SplitObject(srat, split.by = integ_var)
  
  # CCA integration
  features <- SelectIntegrationFeatures(srat_list, nfeatures = 5000)
  srat_list <- PrepSCTIntegration(srat_list, anchor.features = features)
  anchors <- FindIntegrationAnchors(
    srat_list,
    normalization.method = "SCT",
    anchor.features = features,
    reduction = "cca"
  )
  integ_srat <- IntegrateData(anchorset = anchors, normalization.method = "SCT")
  
  # PCA
  integ_srat <- RunPCA(integ_srat, npcs = npcs)
  ggsave(
    "pca_elbow.png",
    plot = ElbowPlot(integ_srat, ndims = npcs),
    path = out_dir,
    width = 8,
    height = 6,
    units = "in"
  )
  
  # UMAP
  integ_srat <- RunUMAP(integ_srat, reduction = "pca", dims = 1:ndims)
  
  # clustering
  integ_srat <- FindNeighbors(integ_srat, dims = 1:ndims) %>% 
    FindClusters(.)
  
  return(integ_srat)
}


#' Use Harmony to integrate a Seurat object.
#' 
#' @param srat Seurat object.
#' @param out_dir Directory to save output plots.
#' @param integ_var Metadata column to integrate on.
#' @param npcs Number of principal components to calculate in `RunPCA`.
#' @param ndims Number of principal components to use for downstream analysis.
#'
#' @return Seurat object after Harmony integration.
#'
integ_harmony <- function(
  srat,
  out_dir,
  integ_var = "orig.ident",
  npcs = 100,
  ndims = 50
) {
  # run sctransform
  srat <- SCTransform(
    srat,
    assay = "RNA",
    variable.features.n = 5000,
    vars.to.regress = "CC.Difference",
    return.only.var.genes = FALSE
  )
  
  # PCA
  srat <- RunPCA(srat, npcs = npcs)
  
  # run harmony
  integ_srat <- RunHarmony(
    object = srat,
    group.by.vars = integ_var,
    plot_convergence = TRUE,
    assay.use = "SCT",
    project.dim = FALSE
  )
  ggsave(
    "pca_elbow.png",
    plot = ElbowPlot(integ_srat, ndims = npcs, reduction = "harmony"),
    path = out_dir,
    width = 8,
    height = 6,
    units = "in"
  )
  
  # UMAP
  integ_srat <- RunUMAP(integ_srat, reduction = "harmony", dims = 1:ndims)
  
  # clustering
  integ_srat <- FindNeighbors(integ_srat, reduction = "harmony", dims = 1:ndims) %>% 
    FindClusters(.)
  
  return(integ_srat)
}


#' Use fastMNN to integrate a Seurat object.
#' 
#' @param srat Seurat object.
#' @param out_dir Directory to save output plots.
#' @param integ_var Metadata column to integrate on.
#' @param npcs Number of principal components to calculate in `RunPCA`.
#' @param ndims Number of principal components to use for downstream analysis.
#'
#' @return Seurat object after fastMNN integration.
#' 
integ_fastmnn <- function(
  srat,
  out_dir,
  integ_var = "orig.ident",
  npcs = 100,
  ndims = 50
) {
  # run sctransform
  srat <- SCTransform(
    srat,
    assay = "RNA",
    variable.features.n = 5000,
    vars.to.regress = "CC.Difference",
    return.only.var.genes = FALSE
  )
  
  # split Seurat into list; `DietSeurat` removes scale.data (otherwise
  # RunFastMNN in the next step gives an error, not sure why)
  srat_list <- DietSeurat(srat) %>% 
    SplitObject(split.by = integ_var)
  
  # run fastMNN
  integ_srat <- SeuratWrappers::RunFastMNN(
    srat_list,
    assay = "SCT",
    features = 5000
  )
  
  # UMAP
  integ_srat <- RunUMAP(integ_srat, reduction = "mnn", dims = 1:ndims)
  
  # clustering
  integ_srat <- FindNeighbors(integ_srat, reduction = "mnn", dims = 1:ndims) %>% 
    FindClusters(.)
  
  return(integ_srat)
}


#' Generate multiple plots with Seurat `DimPlot`.
#'
#' @param srat Seurat object.
#' @param reduction Any dimensional reduction present in the Seurat object.
#' @param group.by Metadata variables to group by.
#' @param label Vector of Booleans, whether or not to add labels to each of the
#'   plots.
#' @param legend Vector of Booleans, whether or not to add a legend to each of
#'   the plots.
#'
#' @return List of `ggplot` objects.
#'
make_plts <- function(
  srat,
  reduction,
  group.by,
  label,
  legend
) {
  plt_list <- purrr::pmap(
    .l = list(group.by, label, legend),
    .f = \(group.by, label, legend) {
      plt <- DimPlot(
        srat,
        reduction = reduction,
        group.by = group.by,
        label = label
      )
      
      if (!legend) {
        plt <- plt + NoLegend()
      }
      
      return(plt)
    }
  )
  
  return(plt_list)
}

# ------------------------------------------------------------------------------
# run the different integration methods

cca_dir <- file.path(args$out_dir, "cca")
harmony_dir <- file.path(args$out_dir, "harmony")
mnn_dir <- file.path(args$out_dir, "mnn")

for (d in c(cca_dir, harmony_dir, mnn_dir)) {
  if (!dir.exists(d)) {
    dir.create(d, recursive = TRUE)
  }
}

# CCA integration
message("\n***RUNNING CCA INTEGRATION***\n")
cca_srat <- integ_cca(srat, out_dir = cca_dir, integ_var = "orig.ident")
saveRDS(cca_srat, file = "mb_cca.rds")

# Harmony integration
message("\n***RUNNING HARMONY INTEGRATION***\n")
harmony_srat <- integ_harmony(srat, out_dir = harmony_dir, integ_var = "orig.ident")
saveRDS(harmony_srat, file = "mb_harmony.rds")

# fastMNN integration
message("\n***RUNNING MNN INTEGRATION***\n")
mnn_srat <- integ_fastmnn(srat, out_dir = mnn_dir, integ_var = "orig.ident")
saveRDS(mnn_srat, file = "mb_mnn.rds")

# ------------------------------------------------------------------------------

# plot cca
purrr::walk(
  .x = c("pca", "umap"),
  .f = \(reduction) {
    make_plts(
      cca_srat,
      reduction = reduction,
      group.by = c("orig.ident", "seurat_clusters", "subtype"),
      label = c(FALSE, TRUE, TRUE),
      legend = c(TRUE, FALSE, FALSE)
    ) %>% 
      wrap_plots() %>% 
      ggsave(
        plot = .,
        filename = paste0(reduction, ".png"),
        path = cca_dir,
        width = 16,
        height = 4,
        units = "in",
        dpi = 600
      )
  }
)

# plot harmony
purrr::walk(
  .x = c("harmony", "umap"),
  .f = \(reduction) {
    make_plts(
      harmony_srat,
      reduction = reduction,
      group.by = c("orig.ident", "seurat_clusters", "subtype"),
      label = c(FALSE, TRUE, TRUE),
      legend = c(TRUE, FALSE, FALSE)
    ) %>% 
      wrap_plots() %>% 
      ggsave(
        plot = .,
        filename = paste0(reduction, ".png"),
        path = harmony_dir,
        width = 16,
        height = 4,
        units = "in",
        dpi = 600
      )
  }
)

# plot fastmnn
purrr::walk(
  .x = c("mnn", "umap"),
  .f = \(reduction) {
    make_plts(
      mnn_srat,
      reduction = reduction,
      group.by = c("orig.ident", "seurat_clusters", "subtype"),
      label = c(FALSE, TRUE, TRUE),
      legend = c(TRUE, FALSE, FALSE)
    ) %>% 
      wrap_plots() %>% 
      ggsave(
        plot = .,
        filename = paste0(reduction, ".png"),
        path = mnn_dir,
        width = 16,
        height = 4,
        units = "in",
        dpi = 600
      )
  }
)

# ------------------------------------------------------------------------------

message("\n***SESSION INFO***\n")
print(sessionInfo())

