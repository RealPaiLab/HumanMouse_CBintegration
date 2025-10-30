# ==============================================================================
# Run Monocle 3 on integrated human UBC subclusters, using CytoTRACE scores to
# set the root.
# ==============================================================================

library(argparse)
library(tidyverse)
library(Seurat)
library(SeuratWrappers)
library(monocle3)

# parse command line arguments
parser <- ArgumentParser()

parser$add_argument(
  # differential gene expression results
  "--srat_file",
  default = NULL,
  required = TRUE
)
parser$add_argument(
  # metadata column containing UBC cluster information
  "--cluster",
  default = NULL,
  required = TRUE
)
parser$add_argument(
  # cytotrace results (raw cytotrace output)
  "--ct_rds",
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
    "--srat_file", "/.mounts/labs/pailab/public/HumanMouseUBC/data/UBC.Harmony.RDS",
    "--cluster", "SCT_snn_res.0.5",
    "--ct_rds", "/.mounts/labs/pailab/private/projects/HumanMouseUBC/integrated_human_ubc/20250312/with_dataset_correction/cytotrace_output.rds",
    "--out_dir", "/.mounts/labs/pailab/private/projects/HumanMouseUBC/integrated_human_ubc/20250313/"
  ))
} else {
  arg_list <- parser$parse_args()
}

message(sprintf("Saving files to %s", arg_list$out_dir))
if (!dir.exists(arg_list$out_dir)) {
  dir.create(arg_list$out_dir, recursive = TRUE)
}

# ------------------------------------------------------------------------------
# functions, set variables

# set default `cell_size` and `cell_stroke` when plotting to get rid of the ugly
# black border that Monocle makes
cell_size <- 0
cell_stroke <- 0.75

source("/.mounts/labs/pailab/private/icheong/CBL_scRNAseq/software/utilities/monocle_helpers.R")

# ------------------------------------------------------------------------------
# load Seurat object and CytoTRACE results

srat <- readRDS(arg_list$srat_file)
srat$dataset_name <- str_remove(srat$dataset_name, "_full_cerebellum_human")
srat$age <- factor(srat$age) %>%
  fct_relevel(\(x) {str_sort(x, numeric = TRUE)})

# use variable features for Monocle
var_features <- VariableFeatures(srat, assay = "integrated")

# load CytoTRACE results
ct_res <- readRDS(arg_list$ct_rds)

# ------------------------------------------------------------------------------
# run Monocle 3 pseudotime

# convert to `cell_data_set` object
# use "integrated" assay (see notes from 2024-06-11)
ubc_cds <- as.cell_data_set(srat, assay = "integrated") %>%
  estimate_size_factors()

# normalize and calculate PCA (`preprocess_cds`)
message("***Running PCA (`preprocess_cds`)***")
ubc_cds <- preprocess_cds(
  ubc_cds,
  method = "PCA", # default is PCA
  num_dim = 100, # default is 50
  norm_method = "none", # "integrated" assay already returns "normalized counts"
  use_genes = var_features
)

# show proportion of variance explained
ggsave(
  "elbow_plot.png",
  plot = plot_pc_variance_explained(ubc_cds),
  path = arg_list$out_dir,
  width = 8,
  height = 8,
  units = "in",
  dpi = 600
)

# keep only top 25 PCs
reducedDims(ubc_cds)[["PCA"]] <- reducedDims(ubc_cds)[["PCA"]][, 1:25]

# dimensionality reduction to UMAP
message("***Running UMAP (`reduce_dimension`)***")
ubc_cds <- reduce_dimension(
  ubc_cds,
  reduction_method = "UMAP",
  preprocess_method = "PCA",
  umap.min_dist = 0.3,
  umap.n_neighbor = 30L
)

# cluster cells with defaults
message("***Grouping cells into clusters (`cluster_cells`)***")
ubc_cds <- cluster_cells(ubc_cds, reduction_method = "UMAP")

# plot clusters, partitions, original UBC clusters, dataset, and age
walk(
  .x = c(
    "cluster",
    "partition",
    arg_list$cluster,
    "dataset_name",
    "age"
  ),
  .f = \(x) {
    alpha <- 1
    title <- x

    # custom settings for different plots
    if (x == "dataset_name") {
      # make points semi-transparent for dataset
      alpha <- 0.25
      title <- "dataset"
    } else if (x == arg_list$cluster) {
      title <- "orig_cluster"
    } else if (x %in% c("cluster", "partition")) {
      title <- paste0("monocle_", x)
    }

    # make plot
    .plt <- plot_cells(
      ubc_cds,
      reduction = "UMAP",
      color_cells_by = x,
      group_cells_by = x,
      label_cell_groups = FALSE,
      show_trajectory_graph = FALSE,
      cell_size = cell_size,
      cell_stroke = cell_stroke,
      alpha = alpha
    )

    .plt <- .plt + 
      labs(title = title)

    ggsave(
      filename = paste0("umap_", title, ".png"),
      path = arg_list$out_dir,
      width = 5,
      height = 4,
      units = "in",
      dpi = 600
    )
  }
)

# learn trajectory
ubc_cds <- learn_graph(ubc_cds, use_partition = FALSE)

# plot cells by trajectory
.plt <- plot_cells(
  ubc_cds,
  color_cells_by = arg_list$cluster,
  label_cell_groups = FALSE,
  cell_size = cell_size,
  cell_stroke = cell_stroke,
  label_principal_points = TRUE
)
ggsave(
  "umap_trajectory.png",
  plot = .plt,
  path = arg_list$out_dir,
  width = 5,
  height = 4,
  units = "in",
  dpi = 600
)

# save the `cell_data_set` before running pseudotime
qs::qsave(
  x = ubc_cds,
  file = file.path(arg_list$out_dir, "pre_pseudotime_cds.qs")
)

# order the cells by pseudotime, setting the root ("time 0") as the cell with
# the highest CytoTRACE score (meaning it's the least differentiated)
starting_cells <- names(ct_res$CytoTRACE[which.max(ct_res$CytoTRACE)])
ubc_cds <- order_cells(
  ubc_cds,
  root_cells = starting_cells
)

# plot cells by pseudotime
.plt <- plot_cells(
  ubc_cds,
  color_cells_by = "pseudotime",
  label_cell_groups = FALSE,
  cell_size = cell_size,
  cell_stroke = cell_stroke
)
ggsave(
  "umap_pseudotime.png",
  plot = .plt,
  path = arg_list$out_dir,
  width = 5,
  height = 4,
  units = "in",
  dpi = 600
)

# save the `cell_data_set` object
qs::qsave(
  x = ubc_cds,
  file = file.path(arg_list$out_dir, "post_pseudotime_cds.qs")
)

# save the pseudotime value of each cell in a separate file for easier
# downstream analysis
save_pseudotime_to_csv(
  cds = ubc_cds,
  file = file.path(arg_list$out_dir, "monocle3_pseudotime.csv")
)

# ------------------------------------------------------------------------------
# comparing Monocle 3 pseudotime with CytoTRACE

# add CytoTRACE scores to `cell_data_set` object
cell_ids <- rownames(colData(ubc_cds))
ubc_cds$cytotrace <- ct_res$CytoTRACE[cell_ids]

# plot CytoTRACE scores on UMAP
.plt <- plot_cells(
  ubc_cds,
  color_cells_by = "cytotrace",
  label_cell_groups = FALSE,
  cell_size = cell_size,
  cell_stroke = cell_stroke
)
ggsave(
  "umap_cytotrace.png",
  plot = .plt,
  path = arg_list$out_dir,
  width = 5,
  height = 4,
  units = "in",
  dpi = 600
)

# get Monocle 3 and CytoTRACE scores, merge into a dataframe
pt_scores <- data.frame(
  cell_id = cell_ids,
  monocle3 = pseudotime(ubc_cds)[cell_ids],
  cytotrace = ct_res$CytoTRACE[cell_ids]
) %>%
  # scale Monocle 3 pseudotime values from 0-1 and plot against CytoTRACE scores
  mutate(monocle3 = scales::rescale(x = monocle3))

# plot Monocle 3 vs CytoTRACE values
.plt <- ggplot(pt_scores, aes(x = monocle3, y = cytotrace)) + 
  geom_point(alpha = 0.5) + 
  theme_light()
ggsave(
  file = "scatter_monocle3_vs_cytotrace.png",
  plot = .plt,
  path = arg_list$out_dir,
  width = 6,
  height = 6,
  units = "in",
  dpi = 600
)

# ------------------------------------------------------------------------------

message("\n***SESSION INFO***\n")
print(sessionInfo())
