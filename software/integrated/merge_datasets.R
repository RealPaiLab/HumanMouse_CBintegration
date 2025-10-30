# ==============================================================================
# Merge huamn and mouse datasets without integrating to visualize the batch
# effect and demonstrate the need for integration.
# ==============================================================================

library(argparse)
library(tidyverse)
library(Seurat)

# parse command line arguments
parser <- ArgumentParser()

parser$add_argument(
  # output directory
  "--out_dir",
  default = NULL,
  required = TRUE
)

if (interactive()) {
  arg_list <- parser$parse_args(c(
    "--out_dir", "/.mounts/labs/pailab/private/icheong/CBL_scRNAseq/results/integrated/20241223"
  ))
} else {
  arg_list <- parser$parse_args()
}

message(sprintf("***Saving files to %s***", arg_list$out_dir))
if (!dir.exists(arg_list$out_dir)) {
  dir.create(arg_list$out_dir, recursive = TRUE)
}

# ------------------------------------------------------------------------------
# functions

source("/.mounts/labs/pailab/private/icheong/CBL_scRNAseq/thesis/utils.R")

my_pals <- get_custom_pals()

# ------------------------------------------------------------------------------
# load Seurat objects

# functions from `utils.R`
srat_qs <- get_srat_paths()
srat_list <- load_srat(
  srat_qs[c("aldinger_human", "sepp_human", "vladoiu_mouse", "sepp_mouse")]
)

# ------------------------------------------------------------------------------
# merge the datasets (following this workflow:
# https://github.com/satijalab/seurat/issues/5738#issuecomment-1169321584)

message("***Merging datasets***")
start <- Sys.time()
srat_mrg <- merge(
  x = srat_list[[1]],
  y = as.vector(srat_list[-1]),
  merge.data = TRUE, # merge data slot
  merge.dr = FALSE # don't merge dimensionality reductions
)
end <- Sys.time()
message(sprintf(
  "***Merging time: %s min***",
  difftime(end, start, units = "mins")
))

# set variable features after merging (see
# https://github.com/satijalab/seurat/issues/6185)
VariableFeatures(srat_mrg) <- SelectIntegrationFeatures(
  object.list = srat_list,
  nfeatures = 5000
)

# correct for any normalization differences between datasets (see
# https://github.com/satijalab/seurat/issues/7407)
message("***Running PrepSCTFindMarkers***")
start <- Sys.time()
srat_mrg <- PrepSCTFindMarkers(srat_mrg)
end <- Sys.time()
message(sprintf(
  "***PrepSCTFindMarkers time: %s min***",
  difftime(end, start, units = "mins")
))

# ------------------------------------------------------------------------------
# run dimensionality reduction

message("***Running PCA***")
srat_mrg <- RunPCA(srat_mrg, npcs = 100)
.plt <- ElbowPlot(srat_mrg, ndims = 100) + theme_classic()
ggsave(
  filename = "elbow_plot.png",
  plot = .plt,
  path = arg_list$out_dir,
  width = 7,
  height = 5,
  units = "in",
  dpi = 600
)

message("***Running UMAP***")
ndims <- 50
srat_mrg <- RunUMAP(srat_mrg, reduction = "pca", dims = 1:ndims)

# save results
message("***Saving merged Seurat object***")
qs::qsave(
  x = srat_mrg,
  file = file.path(arg_list$out_dir, "merged_seurat.qs")
)

# ------------------------------------------------------------------------------
# visualize results

# plot by species
walk(
  .x = c("pca", "umap"),
  .f = \(reduction) {
    .plt <- DimPlot(
      srat_mrg,
      reduction = reduction,
      cells.highlight = WhichCells(srat_mrg, expression = species == "human"),
      sizes.highlight = 0.01,
      label = FALSE,
      raster = FALSE
    ) + 
      labs(title = NULL) + 
      scale_colour_manual(
        labels = c("mouse", "human"),
        values = my_pals$species
      )
    ggsave(
      filename = paste0("species_", reduction, ".png"),
      plot = .plt,
      path = arg_list$out_dir,
      width = 5.5,
      height = 5,
      units = "in",
      dpi = 600
    )
  }
)

# plot by dataset
srat_mrg$dataset_name <- srat_mrg$dataset_name %>%
  str_remove(pattern = "full_cerebellum_")
walk(
  .x = c("pca", "umap"),
  .f = \(reduction) {
    .plt <- DimPlot(
      srat_mrg,
      reduction = reduction,
      group.by = "dataset_name",
      label = FALSE,
      raster = FALSE
    ) + 
      labs(title = NULL)
    ggsave(
      filename = paste0("dataset_", reduction, ".png"),
      plot = .plt,
      path = arg_list$out_dir,
      width = 6,
      height = 5,
      units = "in",
      dpi = 600
    )
  }
)

# ------------------------------------------------------------------------------

message("\n***SESSION INFO***\n")
print(sessionInfo())
