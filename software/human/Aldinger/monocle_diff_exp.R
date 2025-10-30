# ==============================================================================
# Find differentially expressed genes across pseudotime of UBC lineage using
# Monocle 3.
# ==============================================================================

# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(argparse)
library(tidyverse)
library(monocle3)
# library(patchwork)

# parse command line arguments
parser <- ArgumentParser()

parser$add_argument(
  # monocle cds file
  "--monocle_cds",
  default = NULL,
  required = TRUE
)
parser$add_argument(
  # number of cores to use
  "--num_cores",
  default = 16,
  required = FALSE
)
parser$add_argument(
  # output directory
  "--out_dir",
  default = NULL,
  required = TRUE
)

if (interactive()) {
  arg_list <- parser$parse_args(c(
    "--monocle_cds", "/CBL_scRNAseq/results/human/Aldinger/20230719/monocle_cds_ubc_lineage.rds",
    "--out_dir", "/CBL_scRNAseq/results/human/Aldinger/20230809/"
  ))
} else {
  arg_list <- parser$parse_args()
}

message(sprintf("Saving files to %s", getwd()))
if (!dir.exists(arg_list$out_dir)) {
  dir.create(arg_list$out_dir, recursive = TRUE)
}

# import functions
source("/CBL_scRNAseq/software/utilities/monocle_helpers.R")

# ------------------------------------------------------------------------------
# load data

cds <- readRDS(arg_list$monocle_cds)

# ------------------------------------------------------------------------------
# monocle differential expression and gene modules

purrr::walk(
  .x = c("knn", "principal_graph"),
  .f = \(neighbor_graph) {
    
    subdir <- file.path(arg_list$out_dir, neighbor_graph)
    if (!dir.exists(subdir)) dir.create(subdir, recursive = TRUE)
    
    # find differentially expressed genes (wrapper for `graph_test`)
    gt_res <- run_graph_test(
      cds = cds,
      neighbor_graph = neighbor_graph,
      cores = arg_list$num_cores,
      deg_file = file.path(subdir, "deg.csv")
    )
    
    # get names of differentially expressed genes
    deg_names <- base::rownames(filter(gt_res, q_value < 0.05))
    
    # cluster genes into modules (wrapper for `find_gene_modules`)
    modules_df <- get_modules(
      cds = cds[deg_names, ],
      resolution = 10^seq(-1,-6),
      module_file = file.path(subdir, "module_genes.csv")
    )
    
    # subset cell metadata
    cell_annot_df <- tibble(
      cell = base::rownames(colData(cds)),
      integ_clusters = colData(cds)[, "integ_clusters"],
      cell_type = colData(cds)[, "new_cell_type"]
    )
    
    # make heatmap showing gene module expression with all cell annotations
    purrr::walk(
      .x = seq(2, ncol(cell_annot_df)),
      .f = \(x) {
        module_expr <- get_module_expr(
          cds = cds,
          modules_df = modules_df,
          cell_annot_df = cell_annot_df[, c(1, x)],
          out_dir = subdir
        )
        
        hm <- plot_gene_modules(
          module_expr = module_expr,
          cell_annot = colnames(cell_annot_df)[x],
          out_dir = subdir
        )
      }
    )
    
  }
)

# ------------------------------------------------------------------------------

message("\n***SESSION INFO***\n")
print(sessionInfo())
