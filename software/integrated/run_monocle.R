# ==============================================================================
# Perform trajectory inference on the integrated human + mouse dataset using
# Monocle 3.
# ==============================================================================

# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

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

if (interactive()) {
  arg_list <- parser$parse_args(c(
    "--srat_qs", "/.mounts/labs/pailab/private/llau/results/integrated/20240524/25_pc_without_luo/25_pc_rl.qs",
    "--out_dir", "/.mounts/labs/pailab/private/icheong/CBL_scRNAseq/results/integrated/20240611"
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

source("/.mounts/labs/pailab/private/icheong/CBL_scRNAseq/software/utilities/monocle_helpers.R")

# ------------------------------------------------------------------------------
# load Seurat objects

srat <- qs::qread(arg_list$srat_qs)

# remove control cells leaving only RL lineage cells
# oligodendrocytes - cluster 11/15
# endothelial - cluster 13
# microglia - cluster 14
srat <- subset(
  x = srat,
  subset = snn_res.0.4 %in% c(11, 13, 14, 15),
  invert = TRUE
)

# set the resolution as the seurat clusters
srat$seurat_clusters <- srat$snn_res.0.4 %>%
  fct_drop()

# use variable features for Monocle
var_features <- VariableFeatures(srat, assay = "integrated")

# set default `cell_size` and `cell_stroke` when plotting to get rid of the ugly
# black border that Monocle makes
cell_size <- 0
cell_stroke <- 0.6

# also need the UBC subclusters to copy them over
ubc_srat <- qs::qread("/.mounts/labs/pailab/private/llau/results/integrated/20240527/ubc_subset.qs")

# add UBC subclusters to the metadata
md <- merge(
  x = srat[[]],
  y = dplyr::select(ubc_srat[[]], snn_res.0.3) %>%
    dplyr::rename(ubc_subclusters = snn_res.0.3),
  by = "row.names",
  all.x = TRUE
) %>%
  mutate(
    # combine UBC subclusters with other clusters
    final_clusters = case_when(
      !is.na(ubc_subclusters) ~ paste0("UBC_", ubc_subclusters),
      .default = seurat_clusters
    ),
    # set factor levels for the new column
    final_clusters = fct_relevel(
      final_clusters,
      str_sort(names(table(final_clusters)), numeric = TRUE)
    )
  ) %>%
  column_to_rownames(var = "Row.names")

# add cell type annotations to the metadata
cluster_annot <- read_csv("/.mounts/labs/pailab/private/llau/results/integrated/20240618/clusters_cell_types.csv") %>%
  mutate(
    final_clusters = str_remove(Cluster, "^Cluster_"),
    .before = 1
  ) %>%
  dplyr::rename(cell_type_annot = broad_w_ubc_subtypes) %>%
  dplyr::select(final_clusters, cell_type_annot)
md <- merge(
  x = rownames_to_column(md, var = "cell_id"),
  y = cluster_annot,
  by = "final_clusters",
  all.x = TRUE
) %>%
  column_to_rownames(var = "cell_id")

# keep same order as before and replace the metadata
md <- md[rownames(srat[[]]), ]
srat@meta.data <- md

# ------------------------------------------------------------------------------
# load cells of interest (for plotting purposes)

source("/.mounts/labs/pailab/private/icheong/CBL_scRNAseq/software/utilities/load_cell_ids.R")
cell_ids <- load_cell_ids()

# ------------------------------------------------------------------------------
# pseudotime on human cells only

sub_dir <- "human_pseudotime"
message(sprintf(
  "***Running pseudotime on human cells, saving files to %s***",
  file.path(arg_list$out_dir, sub_dir)
))

# convert to `cell_data_set` object
human_cds <- subset(
  x = srat,
  subset = species == "human"
) %>%
  as.cell_data_set(assay = "integrated") %>%
  estimate_size_factors()

# normalize and calculate PCA (`preprocess_cds`)
message("***Running PCA (`preprocess_cds`)***")
human_cds <- preprocess_cds(
  human_cds,
  method = "PCA", # default is PCA
  num_dim = 100, # default is 50
  norm_method = "none", # "integrated" assay already returns "normalized counts"
  use_genes = var_features
)

# show proportion of variance explained
ggsave(
  "elbow_plot.png",
  plot = plot_pc_variance_explained(human_cds),
  path = file.path(arg_list$out_dir, sub_dir),
  width = 8,
  height = 8,
  units = "in",
  dpi = 600
)

# keep only top 25 PCs (see notebook)
reducedDims(human_cds)[["PCA"]] <- reducedDims(human_cds)[["PCA"]][, 1:25]

# dimensionality reduction to UMAP
message("***Running UMAP (`reduce_dimension`)***")
human_cds <- reduce_dimension(
  human_cds,
  reduction_method = "UMAP",
  preprocess_method = "PCA",
  umap.min_dist = 0.3,
  umap.n_neighbor = 30L
)

# cluster cells with defaults
message("***Grouping cells into clusters (`cluster_cells`)***")
human_cds <- cluster_cells(human_cds, reduction_method = "UMAP")

# plot clusters, partitions, original Seurat clusters, and cell type
walk(
  .x = c("cluster", "partition", "seurat_clusters", "cell_type_annot"),
  .f = \(x) {
    # make plot
    .plt <- plot_cells(
      human_cds,
      reduction = "UMAP",
      color_cells_by = x,
      group_cells_by = x,
      label_cell_groups = FALSE,
      show_trajectory_graph = FALSE,
      cell_size = cell_size,
      cell_stroke = cell_stroke
    )

    # add title to plot
    if (x %in% c("cluster", "partition")) {x <- paste0("monocle_", x)}
    .plt <- .plt + 
      labs(title = x)

    ggsave(
      filename = paste0("umap_", x, ".png"),
      plot = .plt,
      path = file.path(arg_list$out_dir, sub_dir),
      width = 7,
      height = 6,
      units = "in",
      dpi = 600
    )
  }
)

# highlight some cells of interest on the UMAP (RL-VZ/SVZ, UBCs from previous
# integration)
.plt1 <- plot_highlight(
  cds = human_cds,
  cells = cell_ids[c("rl_vz", "rl_svz")],
  highlight_cols = RColorBrewer::brewer.pal(n = 5, name = "Set1")[1:2]
) + 
  labs(title = "RL-VZ and -SVZ cells (Hendrikse)")
.plt2 <- plot_highlight(
  cds = human_cds,
  cells = cell_ids[paste0("cluster_", c(7, 19, 20))],
  highlight_cols = RColorBrewer::brewer.pal(n = 5, name = "Set1")[3:5]
) + 
  labs(title = "UBC clusters from initial integration")
.plt <- (.plt1 + .plt2) & 
  guides(colour = guide_legend(nrow = 1, override.aes = list(size = 2))) &
  theme(legend.position = "bottom")
ggsave(
  filename = "umap_orig_cells.png",
  plot = .plt,
  path = file.path(arg_list$out_dir, sub_dir),
  width = 10,
  height = 6,
  units = "in",
  dpi = 600
)

# learn trajectory
human_cds <- learn_graph(human_cds, use_partition = FALSE)

# plot cells by trajectory
.plt <- plot_cells(
  human_cds,
  color_cells_by = "cell_type_annot",
  label_cell_groups = FALSE,
  cell_size = cell_size,
  cell_stroke = cell_stroke,
  label_principal_points = TRUE
)
ggsave(
  "umap_trajectory.png",
  plot = .plt,
  path = file.path(arg_list$out_dir, sub_dir),
  width = 7,
  height = 6,
  units = "in",
  dpi = 600
)

# save the `cell_data_set` before running pseudotime
qs::qsave(
  x = human_cds,
  file = file.path(arg_list$out_dir, sub_dir, "pre_pseudotime_cds.qs")
)

# order the cells by pseudotime, setting the root ("time 0") as the RL cells
# (Seurat cluster 10)
starting_cells <- WhichCells(
  srat,
  expression = species == "human" & seurat_clusters == 10
)
human_cds <- order_cells(
  human_cds,
  root_cells = starting_cells
)

# plot cells by pseudotime
.plt <- plot_cells(
  human_cds,
  color_cells_by = "pseudotime",
  label_cell_groups = FALSE,
  cell_size = cell_size,
  cell_stroke = cell_stroke
)
ggsave(
  "umap_pseudotime.png",
  plot = .plt,
  path = file.path(arg_list$out_dir, sub_dir),
  width = 7,
  height = 6,
  units = "in",
  dpi = 600
)

# save the `cell_data_set` object after running pseudotime
qs::qsave(
  x = human_cds,
  file = file.path(arg_list$out_dir, sub_dir, "post_pseudotime_cds.qs")
)

# save the pseudotime value of each cell in a separate file for easier
# downstream analysis
save_pseudotime_to_csv(
  cds = human_cds,
  file = file.path(arg_list$out_dir, sub_dir, "pseudotime_values.csv")
)

# ------------------------------------------------------------------------------
# pseudotime on mouse cells only

sub_dir <- "mouse_pseudotime"
message(sprintf(
  "***Running pseudotime on mouse cells, saving files to %s***",
  file.path(arg_list$out_dir, sub_dir)
))

# convert to `cell_data_set` object
mouse_cds <- subset(
  x = srat,
  subset = species == "mouse"
) %>%
  as.cell_data_set(assay = "integrated") %>%
  estimate_size_factors()

# normalize and calculate PCA (`preprocess_cds`)
message("***Running PCA (`preprocess_cds`)***")
mouse_cds <- preprocess_cds(
  mouse_cds,
  method = "PCA", # default is PCA
  num_dim = 100, # default is 50
  norm_method = "none", # "integrated" assay already returns "normalized counts"
  use_genes = var_features
)

# show proportion of variance explained
ggsave(
  "elbow_plot.png",
  plot = plot_pc_variance_explained(mouse_cds),
  path = file.path(arg_list$out_dir, sub_dir),
  width = 8,
  height = 8,
  units = "in",
  dpi = 600
)

# keep only top 25 PCs (see notebook)
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
  .x = c("cluster", "partition", "seurat_clusters", "cell_type_annot"),
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
      path = file.path(arg_list$out_dir, sub_dir),
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
  color_cells_by = "cell_type_annot",
  label_cell_groups = FALSE,
  cell_size = cell_size,
  cell_stroke = cell_stroke,
  label_principal_points = TRUE
)
ggsave(
  "umap_trajectory.png",
  plot = .plt,
  path = file.path(arg_list$out_dir, sub_dir),
  width = 7,
  height = 6,
  units = "in",
  dpi = 600
)

# save the `cell_data_set` before running pseudotime
qs::qsave(
  x = mouse_cds,
  file = file.path(arg_list$out_dir, sub_dir, "pre_pseudotime_cds.qs")
)

# order the cells by pseudotime, setting the root ("time 0") as the RL cells
# (Seurat cluster 12)
starting_cells <- WhichCells(
  srat,
  expression = species == "mouse" & seurat_clusters == 12
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
  cell_stroke = cell_stroke
)
ggsave(
  "umap_pseudotime.png",
  plot = .plt,
  path = file.path(arg_list$out_dir, sub_dir),
  width = 7,
  height = 6,
  units = "in",
  dpi = 600
)

# save the `cell_data_set` object
qs::qsave(
  x = mouse_cds,
  file = file.path(arg_list$out_dir, sub_dir, "post_pseudotime_cds.qs")
)

# save the pseudotime value of each cell in a separate file for easier
# downstream analysis
save_pseudotime_to_csv(
  cds = mouse_cds,
  file = file.path(arg_list$out_dir, sub_dir, "pseudotime_values.csv")
)

# ------------------------------------------------------------------------------
# pseudotime on integrated human and mouse cells

sub_dir <- "integ_pseudotime"
message(sprintf(
  "***Running pseudotime on integrated human and mouse cells, saving files to %s***",
  file.path(arg_list$out_dir, sub_dir)
))

# convert to `cell_data_set` object; also remove cluster 12 as this appears to
# be a mouse-only RL or pre-RL cluster from an earlier timepoint
integ_cds <- subset(
  x = srat,
  subset = snn_res.0.4 %in% c(12),
  invert = TRUE
) %>%
  as.cell_data_set(assay = "integrated") %>%
  estimate_size_factors()

# normalize and calculate PCA (`preprocess_cds`)
message("***Running PCA (`preprocess_cds`)***")
integ_cds <- preprocess_cds(
  integ_cds,
  method = "PCA", # default is PCA
  num_dim = 100, # default is 50
  norm_method = "none", # "integrated" assay already returns "normalized counts"
  use_genes = var_features
)

# show proportion of variance explained
ggsave(
  "elbow_plot.png",
  plot = plot_pc_variance_explained(integ_cds),
  path = file.path(arg_list$out_dir, sub_dir),
  width = 8,
  height = 8,
  units = "in",
  dpi = 600
)

# keep only top 25 PCs (see notebook)
reducedDims(integ_cds)[["PCA"]] <- reducedDims(integ_cds)[["PCA"]][, 1:25]

# dimensionality reduction to UMAP
message("***Running UMAP (`reduce_dimension`)***")
integ_cds <- reduce_dimension(
  integ_cds,
  reduction_method = "UMAP",
  preprocess_method = "PCA",
  umap.min_dist = 0.3,
  umap.n_neighbor = 30L
)

# cluster cells with defaults
message("***Grouping cells into clusters (`cluster_cells`)***")
integ_cds <- cluster_cells(integ_cds, reduction_method = "UMAP")

# plot clusters, partitions, original Seurat clusters, cell type, and species
walk(
  .x = c(
    "cluster",
    "partition",
    "seurat_clusters",
    "cell_type_annot",
    "species"
  ),
  .f = \(x) {
    # make points semi-transparent for species
    alpha <- 1
    if (x == "species") {
      alpha <- 0.25
    }

    # make plot
    .plt <- plot_cells(
      integ_cds,
      reduction = "UMAP",
      color_cells_by = x,
      group_cells_by = x,
      label_cell_groups = FALSE,
      show_trajectory_graph = FALSE,
      cell_size = cell_size,
      cell_stroke = cell_stroke,
      alpha = alpha
    )

    # add title
    if (x %in% c("cluster", "partition")) {
      x <- paste0("monocle_", x)
    }
    .plt <- .plt + 
      labs(title = x)

    ggsave(
      filename = paste0("umap_", x, ".png"),
      path = file.path(arg_list$out_dir, sub_dir),
      width = 7,
      height = 6,
      units = "in",
      dpi = 600
    )
  }
)

# learn trajectory
integ_cds <- learn_graph(integ_cds, use_partition = FALSE)

# plot cells by trajectory
.plt <- plot_cells(
  integ_cds,
  color_cells_by = "cell_type_annot",
  label_cell_groups = FALSE,
  cell_size = cell_size,
  cell_stroke = cell_stroke,
  label_principal_points = TRUE
)
ggsave(
  "umap_trajectory.png",
  plot = .plt,
  path = file.path(arg_list$out_dir, sub_dir),
  width = 7,
  height = 6,
  units = "in",
  dpi = 600
)

# save the `cell_data_set` before running pseudotime
qs::qsave(
  x = integ_cds,
  file = file.path(arg_list$out_dir, sub_dir, "pre_pseudotime_cds.qs")
)

# order the cells by pseudotime, setting the root ("time 0") as the RL cells
# (Seurat cluster 10)
starting_cells <- WhichCells(
  srat,
  expression = species == "mouse" & seurat_clusters == 10
)
integ_cds <- order_cells(
  integ_cds,
  root_cells = starting_cells
)

# plot cells by pseudotime
.plt <- plot_cells(
  integ_cds,
  color_cells_by = "pseudotime",
  label_cell_groups = FALSE,
  cell_size = cell_size,
  cell_stroke = cell_stroke
)
ggsave(
  "umap_pseudotime.png",
  plot = .plt,
  path = file.path(arg_list$out_dir, sub_dir),
  width = 7,
  height = 6,
  units = "in",
  dpi = 600
)

# save the `cell_data_set` object
qs::qsave(
  x = integ_cds,
  file = file.path(arg_list$out_dir, sub_dir, "post_pseudotime_cds.qs")
)

# save the pseudotime value of each cell in a separate file for easier
# downstream analysis
save_pseudotime_to_csv(
  cds = integ_cds,
  file = file.path(arg_list$out_dir, sub_dir, "pseudotime_values.csv")
)

# ------------------------------------------------------------------------------

message("\n***SESSION INFO***\n")
print(sessionInfo())
