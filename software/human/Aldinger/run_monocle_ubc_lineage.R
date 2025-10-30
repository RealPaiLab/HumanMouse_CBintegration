# ==============================================================================
# Subset the human RL to UBC lineage of cells and re-run Monocle 3
# ==============================================================================

# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(argparse)
library(tidyverse)
library(monocle3)
library(patchwork)

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
    "--out_dir", "/CBL_scRNAseq/results/human/Aldinger/20230719/"
  ))
} else {
  arg_list <- parser$parse_args()
}

message(sprintf("Saving files to %s", getwd()))
if (!dir.exists(arg_list$out_dir)) {
  dir.create(arg_list$out_dir, recursive = TRUE)
}

# ------------------------------------------------------------------------------

# load Monocle CDS object
cds <- readRDS("/CBL_scRNAseq/results/human/Aldinger/20230718/monocle_cds.rds")

# plot the principal nodes of the full RL
.plt <- plot_cells(
  cds,
  reduction_method = "UMAP",
  color_cells_by = "new_cell_type",
  label_cell_groups = FALSE,
  label_principal_points = TRUE
)
ggsave(
  "full_rl_principal_points.png",
  plot = .plt,
  path = arg_list$out_dir,
  width = 6,
  height = 4,
  units = "in",
  dpi = 600
)

# subset cells based on plot
cds_ubc <- choose_graph_segments(
  cds = cds,
  starting_pr_node = "Y_9",
  ending_pr_nodes = "Y_104"
)

# ------------------------------------------------------------------------------
# redo Monocle 3 PCA and dimensionality reduction

# normalize and PCA
message("Running PCA (`preprocess_cds`)")
cds_ubc <- preprocess_cds(
  cds_ubc,
  method = "PCA",
  num_dim = 100,
  norm_method = "log"
)

# dimensionality reduction with UMAP
message("Running UMAP (`reduce_dimension`)")
cds_ubc <- reduce_dimension(
  cds_ubc,
  reduction_method = "UMAP",
  preprocess_method = "PCA"
)

# ------------------------------------------------------------------------------
# clustering and pseudotime

cds_ubc <- cds_ubc %>% 
  cluster_cells(reduction_method = "UMAP") %>% 
  learn_graph(use_partition = TRUE, close_loop = FALSE) %>% 
  order_cells(
    # set root (time = 0) to RL-VZ cells
    root_cells = rownames(subset(
      x = colData(cds_ubc),
      subset = new_cell_type %in% c("RL-VZ")
    ))
  )



# https://stackoverflow.com/questions/11880906/pass-subset-argument-through-a-function-to-subset <-- continue here!!



# ------------------------------------------------------------------------------
# plotting

message("Saving plots")

# show proportion of variance explained
write_lines(
  x = cds_ubc@reduce_dim_aux$PCA$model$prop_var_expl,
  file = file.path(arg_list$out_dir, "pc_variance_explained.txt")
)
ggsave(
  "pc_variance_explained.png",
  plot = plot_pc_variance_explained(cds),
  path = arg_list$out_dir,
  width = 8,
  height = 4,
  units = "in",
  dpi = 600
)

# plot cells with different colours
.plt <- map(
  .x = c("pseudotime", "new_cell_type", "integ_clusters"),
  .f = \(x) {
    plot_cells(
       cds_ubc,
       color_cells_by = x,
       cell_size = 0.75,
       label_cell_groups = FALSE
    )
  }
) %>% 
  wrap_plots(ncol = 3)
.plt[[3]] <- .plt[[3]] + 
  scale_colour_manual(values = Seurat::DiscretePalette(
    n = length(levels(colData(cds_ubc)$integ_clusters)),
    palette = "polychrome"
  ))
ggsave(
  "trajectory_umaps.png",
  plot = .plt,
  path = arg_list$out_dir,
  width = 16,
  height = 4,
  units = "in"
)

# ------------------------------------------------------------------------------

message(sprintf("Saving %s", file.path(arg_list$out_dir, "monocle_cds_ubc_lineage.rds")))
saveRDS(cds_ubc, file = file.path(arg_list$out_dir, "monocle_cds_ubc_lineage.rds"))

message("\n***SESSION INFO***\n")
print(sessionInfo())
