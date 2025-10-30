# ==============================================================================
# Perform trajectory inference on the mouse datasets using Monocle 3.
# ==============================================================================

library(argparse)
library(tidyverse)
library(patchwork)
library(Seurat)
library(SeuratWrappers)
library(monocle3)

# parse command line arguments
parser <- ArgumentParser()

parser$add_argument(
  # integrated Seurat object
  "--srat_qs",
  default = NULL,
  required = TRUE
)
parser$add_argument(
  # output directory
  "--out_dir",
  default = NULL,
  required = TRUE
)
parser$add_argument(
  # root cells for pseudotime (pass in the metadata column name and a value in the column)
  "--starting_cells",
  nargs = "*",
  required = TRUE
)
parser$add_argument(
  # cells to remove (pass in the metadata column name and a value in the column)
  "--filter_cells",
  nargs = "*",
  required = FALSE
)

if (interactive()) {
  arg_list <- parser$parse_args(c(
    "--srat_qs", "/.mounts/labs/pailab/private/llau/data/Sepp_2024/Sepp_RL_mouse_20240415_stand.qs",
    "--out_dir", "/.mounts/labs/pailab/private/icheong/CBL_scRNAseq/results/mouse/Sepp/20240628",
    "--starting_cells", "cell_type", "GC/UBC_Sepp"#,
    # "--filter_cells", "cell_type", "Upper rhombic lip progenitors_Vladoiu"
  ))
} else {
  arg_list <- parser$parse_args()
}

message(sprintf("***Saving files to %s***", arg_list$out_dir))
if (!dir.exists(arg_list$out_dir)) {
  dir.create(arg_list$out_dir, recursive = TRUE)
}

# ------------------------------------------------------------------------------
# load Seurat object

srat <- qs::qread(arg_list$srat_qs)

# keep RL lineage cells, remove control cells (oligodendrocytes, endothelial,
# microglia)
message("***Removing oligodendrocytes, endothelial cells, and microglia***")
srat <- subset(
  x = srat,
  subset = common_cell_name %in% c("oligodendrocyte/OPC", "endothelial", "microglia"),
  invert = TRUE
)

# remove any other additional cells passed from the command line
if (!is.null(arg_list$filter_cells)) {
  if (length(arg_list$filter_cells) < 2) {
    warning("Both the Seurat metadata column name and cell type/cluster annotation must be passed in via the command line; no cells will be filtered out.")
  } else if (!arg_list$filter_cells[1] %in% colnames(srat[[]])) {
    warning(sprintf(
      "`%s` could not be found in the Seurat metadata; no cells will be filtered out.",
      arg_list$filter_cells[1]
    ))
  } else {
    col_to_filter <- arg_list$filter_cells[1]
    cells_to_filter <- rownames(
      srat[[]][pull(srat[[]], col_to_filter) %in% arg_list$filter_cells[-1], ]
    )
    srat <- subset(
      x = srat,
      cells = cells_to_filter,
      invert = TRUE
    )
  }
}

# use variable features for Monocle
var_features <- VariableFeatures(srat, assay = "SCT")

# set default `cell_size` and `cell_stroke` when plotting to get rid of the ugly
# black border that Monocle makes
cell_size <- 0
cell_stroke <- 0.6

# clean up metadata
srat@meta.data <- mutate(
  srat@meta.data,
  orig_cell_type = str_remove(string = cell_type, pattern = "_([:alpha:])*$")
)

# ------------------------------------------------------------------------------
# running pseudotime

# convert to `cell_data_set` object
mouse_cds <- as.cell_data_set(x = srat, assay = "SCT")

# normalize and calculate PCA
message("***Running PCA (`preprocess_cds`)***")
mouse_cds <- preprocess_cds(
  mouse_cds,
  method = "PCA", # default is PCA
  num_dim = 100, # default is 50
  norm_method = "none", # "SCT" assay already returns normalized counts
  use_genes = var_features
)

# show proportion of variance explained
ggsave(
  "elbow_plot.png",
  plot = plot_pc_variance_explained(mouse_cds),
  path = arg_list$out_dir,
  width = 8,
  height = 8,
  units = "in",
  dpi = 600
)

# keep only top 25 PCs
reducedDims(mouse_cds)[["PCA"]] <- reducedDims(mouse_cds)[["PCA"]][, 1:25]

# dimensionality reduction to UMAP
message("***Running UMAP (`reduce_dimension`)***")
mouse_cds <- reduce_dimension(
  mouse_cds,
  reduction_method = "UMAP",
  preprocess_method = "PCA",
  umap.min_dist = 0.3,
  umap.n_neighbor = 30L
)

# cluster cells with defaults
message("***Grouping cells into clusters (`cluster_cells`)***")
mouse_cds <- cluster_cells(mouse_cds, reduction_method = "UMAP")

# plot clusters, partitions, original Seurat clusters, and cell type
walk(
  .x = c("cluster", "partition", "orig_cell_type", "common_cell_name"),
  .f = \(x) {
    # make plot
    .plt <- plot_cells(
      mouse_cds,
      reduction = "UMAP",
      color_cells_by = x,
      group_cells_by = x,
      label_cell_groups = FALSE,
      show_trajectory_graph = FALSE,
      cell_size = cell_size,
      cell_stroke = cell_stroke
    )

    # add title
    if (x %in% c("cluster", "partition")) {x <- paste0("monocle_", x)}
    .plt <- .plt + 
      labs(title = x)

    ggsave(
      filename = paste0("umap_", x, ".png"),
      path = arg_list$out_dir,
      width = 7,
      height = 6,
      units = "in",
      dpi = 600
    )
  }
)

# learn trajectory
mouse_cds <- learn_graph(mouse_cds, use_partition = FALSE)

# plot cells by trajectory
.plt <- plot_cells(
  mouse_cds,
  label_cell_groups = FALSE,
  cell_size = cell_size,
  cell_stroke = cell_stroke,
  label_principal_points = TRUE
)
ggsave(
  "umap_trajectory.png",
  plot = .plt,
  path = arg_list$out_dir,
  width = 7,
  height = 6,
  units = "in",
  dpi = 600
)

# save the `cell_data_set` before running pseudotime
message("***Saving `pre_pseudotime_cds.qs`***")
qs::qsave(
  x = mouse_cds,
  file = file.path(arg_list$out_dir, "pre_pseudotime_cds.qs")
)

# order the cells by pseudotime
starting_cells <- rownames(
  srat[[]][pull(srat[[]], arg_list$starting_cells[1]) %in% arg_list$starting_cells[-1], ]
)
mouse_cds <- order_cells(
  mouse_cds,
  root_cells = starting_cells
)

# plot cells by pseudotime
.plt <- plot_cells(
  mouse_cds,
  color_cells_by = "pseudotime",
  label_cell_groups = FALSE,
  cell_size = cell_size,
  cell_stroke = cell_stroke,
  label_principal_points = TRUE
)
ggsave(
  "umap_pseudotime.png",
  plot = .plt,
  path = arg_list$out_dir,
  width = 7,
  height = 6,
  units = "in",
  dpi = 600
)

# save the `cell_data_set` object
message("***Saving `post_pseudotime_cds.qs`***")
qs::qsave(
  x = mouse_cds,
  file = file.path(arg_list$out_dir, "post_pseudotime_cds.qs")
)

# ------------------------------------------------------------------------------

message("\n***SESSION INFO***\n")
print(sessionInfo())
