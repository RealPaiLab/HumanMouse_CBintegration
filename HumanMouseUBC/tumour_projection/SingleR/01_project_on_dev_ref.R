# ==============================================================================
# Annotate (project) tumour cells with SingleR using the developing cerebellum
# as a reference.
# ==============================================================================

# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(argparse)
library(tidyverse)
library(patchwork)
library(Seurat)
library(SingleR)

# parse command line arguments
parser <- ArgumentParser()

parser$add_argument(
  # two arguments: Seurat tumour file and column with cluster names
  "--mb_srat",
  nargs = 2,
  default = NULL,
  required = TRUE
)
parser$add_argument(
  # two arguments: Seurat file with human RL lineage cells and column with
  # clusters/cell types
  "--rl_srat",
  nargs = 2,
  default = NULL,
  required = TRUE
)
parser$add_argument(
  # two arguments: Seurat file with human UBC cells and column with
  # clusters/cell types
  "--ubc_srat",
  nargs = 2,
  default = NULL,
  required = TRUE
)
parser$add_argument(
  # number of threads to use for SingleR
  "--num_threads",
  default = 1,
  type = "integer"
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
    "--mb_srat", "/.mounts/labs/pailab/private/icheong/CBL_scRNAseq/results/tumour/Vladoiu/20250211/mb_mnn.qs", "SCT_snn_res.0.25",
    "--rl_srat", "/.mounts/labs/pailab/private/projects/HumanMouseUBC/results/integrated_HsFullCB/20241031/cca/assignClusterIdentity/cca_RLlineage_only_241107.qs", "curtype",
    "--ubc_srat", "/.mounts/labs/pailab/private/projects/HumanMouseUBC/results/integrated_HsFullCB/20241031/cca/RLlineage_only/241107/UBC_withClusterAssignments_241107.qs", "seurat_clusters",
    "--num_threads", "4",
    "--out_dir", "/.mounts/labs/pailab/private/projects/HumanMouseUBC/tumour_projection/SingleR/20250218/"
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

source("/.mounts/labs/pailab/private/icheong/CBL_scRNAseq/software/utilities/cluster_barplot.R")

#' Extract normalized gene expression from Seurat object.
#'
#' @param srat Seurat object.
#' @param rerun_norm Rerun `NormalizeData`? Defaults to `TRUE`.
#'
#' @return Matrix of normalized counts.
#'
get_norm_expr <- function(srat, rerun_norm = TRUE) {
  if (rerun_norm) {
    srat <- NormalizeData(srat, assay = "RNA")
  }
  norm_mat <- GetAssayData(srat, slot = "data", assay = "RNA")
  return(norm_mat)
}

#' Read in Seurat objects from an `rds` or `qs` file.
#'
#' @param path Vector of path(s) to the Seurat objects to read in.
#'
#' @return List of Seurat objects.
#'
load_srat <- function(
  path
) {
  srat <- map(
    .x = path,
    .f = \(pth) {
      message(sprintf("***Reading in %s***", pth))
      if (grepl(pattern = "\\.qs$", x = pth, ignore.case = TRUE)) {
        qs::qread(pth)
      } else if (grepl(pattern = "\\.rds$", x = pth, ignore.case = TRUE)) {
        readRDS(pth)
      }
    }
  )
  if (length(srat) == 1) {
    srat <- pluck(srat, 1)
  }
  return(srat)
}

#' Make and save plots for the SingleR prediction results.
#'
#' @param srat Seurat object of MB tumour.
#' @param srat_pred_col Metadata column in `srat` where predictions
#'   (projections) are stored.
#' @param output Path to save plots to.
#' @param prediction_results Results from running `SingleR::SingleR()` (a
#'   DataFrame).
#' @param annotation_col A data frame passed of additional annotations to add to
#'   the heatmap (e.g., MB tumour subtype). Passed to
#'   `SingleR::plotScoreHeatmap`.
#'
make_plots <- function(
  srat,
  srat_pred_col,
  path,
  prediction_results,
  annotation_col = NULL
) {
  # prediction heatmap
  .plt <- plotScoreHeatmap(
    prediction_results,
    annotation_col = annotation_col,
    cluster_cols = TRUE,
    silent = TRUE
  )
  message("***Saving `score_heatmap.png`***")
  ggsave(
    "score_heatmap.png",
    plot = .plt,
    path = path,
    width = 20,
    height = 8,
    units = "in",
    dpi = 600
  )

  # for distribution plots
  num_cell_types <- length(table(prediction_results$pruned.labels))
  nc <- ceiling(sqrt(num_cell_types)) # number of columns
  wd <- nc * 2 + 1 # total width of plot
  ht <- ceiling(num_cell_types / nc) * 4 # total height of plot

  # distribution of scores
  .plt <- plotScoreDistribution(prediction_results, ncol = nc)
  message("***Saving `score_distr.png`***")
  ggsave(
    "score_distr.png",
    plot = .plt,
    path = path,
    width = wd,
    height = ht,
    units = "in",
    dpi = 600
  )

  # distribution of delta scores from median
  .plt <- plotDeltaDistribution(prediction_results, ncol = nc)
  message("***Saving `delta_distr.png`***")
  ggsave(
    "delta_distr.png",
    plot = .plt,
    path = path,
    width = wd,
    height = ht,
    units = "in",
    dpi = 600
  )


  # barplots
  .plt1 <- cluster_barplot(
    srat,
    split.by = srat_pred_col,
    group.by = "subtype",
    position = "fill"
  )

  .plt2 <- cluster_barplot(
    srat,
    split.by = srat_pred_col,
    group.by = "orig.ident",
    position = "fill"
  ) + 
    theme(axis.text.x = element_text(angle = 30, hjust = 1))

  .plt3 <- cluster_barplot(
    srat,
    split.by = arg_list$mb_srat[2],
    group.by = "subtype",
    position = "fill"
  )

  .plt4 <- cluster_barplot(
    srat,
    split.by = arg_list$mb_srat[2],
    group.by = "orig.ident",
    position = "fill"
  ) + 
    theme(axis.text.x = element_text(angle = 30, hjust = 1))

  .plt5 <- cluster_barplot(
    srat,
    split.by = srat_pred_col,
    group.by = arg_list$mb_srat[2],
    position = "fill",
    width = 0.5
  )

  layout <- c(
    area(t = 1, l = 1, r = 1),
    area(t = 1, l = 2, r = 3),
    area(t = 2, l = 1, r = 3)
  )
  .plt <- (
    .plt3 + .plt4 + 
      plot_layout(widths = c(1, 2), guides = "collect") & 
      guides(fill = guide_legend(ncol = 2))
      # scale_fill_manual(values = DiscretePalette(
      #   length(unique(mb_srat$seurat_clusters)),
      #   "stepped"
      # ))
  ) / 
    (.plt1 + .plt2 + .plt5 + plot_layout(design = layout, guides = "collect") & 
      scale_fill_manual(
        values = DiscretePalette(
          n = length(table(srat[[]][srat_pred_col], useNA = "ifany")),
          palette = "polychrome"
        )
      )
    ) + 
    plot_layout(heights = c(1, 2))
  message("***Saving `annot_bar.png`***")
  ggsave(
    "annot_bar.png",
    plot = .plt,
    path = path,
    width = 8,
    height = 8,
    units = "in",
    dpi = 600
  )


  # UMAPs
  .plt <- purrr::map(
    .x = c("orig.ident", "subtype", srat_pred_col),
    .f = \(group.by) {
      .plt <- DimPlot(
        srat,
        reduction = "umap",
        group.by = group.by,
        label = FALSE
      )
      return(.plt)
    }
  )
  message("***Saving `umaps.png`***")
  ggsave(
    "umaps.png",
    plot = wrap_plots(.plt),
    path = path,
    width = 16,
    height = 4,
    units = "in",
    dpi = 600
  )
}

# ------------------------------------------------------------------------------
# load Seurat objects and extract normalized counts

# `get_norm_expr` re-runs normalization (SingleR needs in normalized raw counts,
# don't use SCTransform) and subsets the normalized expression matrices
# see also https://github.com/SingleR-inc/SingleR/issues/98

# tumour
mb_srat <- load_srat(arg_list$mb_srat[1])
mb_mat <- get_norm_expr(mb_srat, rerun_norm = TRUE)

# RL lineage
rl_srat <- load_srat(arg_list$rl_srat[1])
rl_mat <- get_norm_expr(rl_srat, rerun_norm = TRUE)

# UBC clusters
ubc_srat <- load_srat(arg_list$ubc_srat[1])
ubc_mat <- get_norm_expr(ubc_srat, rerun_norm = TRUE)

# ------------------------------------------------------------------------------
# SingleR with RL lineage as reference

message("***Projecting MB tumour cells on RL lineage cells (+ other controls)***")

subdir <- file.path(arg_list$out_dir, "full_RL_lineage")
if (!dir.exists(subdir)) {
  dir.create(subdir, recursive = TRUE)
}

# get reference annotations (make sure it's in the same order as expr matrix)
rl_labels <- rl_srat[[]][colnames(rl_mat), arg_list$rl_srat[2]]

# note de.method is "wilcox" cuz the reference is a single cell dataset
# (see https://bioconductor.org/packages/release/bioc/vignettes/SingleR/inst/doc/SingleR.html)
annot_rl <- SingleR(
  test = mb_mat,
  ref = rl_mat,
  labels = rl_labels,
  de.method = "wilcox",
  num.threads = arg_list$num_threads
)
message(sprintf(
  "***%s labels were pruned***",
  sum(is.na(annot_rl$pruned.labels))
))

saveRDS(annot_rl, file = file.path(subdir, "prediction_results.rds"))

# add labels back to metadata
mb_srat$full_rl_pred <- annot_rl[rownames(mb_srat[[]]), "pruned.labels"]

# save plots
make_plots(
  srat = mb_srat,
  srat_pred_col = "full_rl_pred",
  path = subdir,
  prediction_results = annot_rl,
  annotation_col = mb_srat[[]]["subtype"]
)

# ------------------------------------------------------------------------------
# SingleR with UBCs as reference; take those that were classified as UBC in the
# previous step and re-classify with UBC subclusters

message("***Projecting MB tumour cells classified as UBCs on UBC subclusters***")

subdir <- file.path(arg_list$out_dir, "full_RL_then_UBC")
if (!dir.exists(subdir)) {
  dir.create(subdir, recursive = TRUE)
}

# get reference annotations
ubc_labels <- ubc_srat[[]][colnames(ubc_mat), arg_list$ubc_srat[2]]
levels(ubc_labels) <- paste0("UBC_", levels(ubc_labels))

# subset tumour cells that were classified as UBCs
subset_cell_ids <- rownames(annot_rl)[annot_rl$pruned.labels == "UBC"]
mb_mat_ubc <- mb_mat[, subset_cell_ids]

# run SingleR
annot_ubc <- SingleR(
  test = mb_mat_ubc,
  ref = ubc_mat,
  labels = ubc_labels,
  de.method = "wilcox",
  num.threads = arg_list$num_threads
)
message(sprintf(
  "***%s labels were pruned***",
  sum(is.na(annot_ubc$pruned.labels))
))

saveRDS(annot_ubc, file = file.path(subdir, "prediction_results.rds"))

# add labels back to metadata
mb_srat$rl_then_ubc_pred <- annot_ubc$pruned.labels[match(rownames(mb_srat[[]]), rownames(annot_ubc), nomatch = NA)]

# save plots
make_plots(
  srat = subset(mb_srat, cells = subset_cell_ids),
  srat_pred_col = "rl_then_ubc_pred",
  path = subdir,
  prediction_results = annot_ubc,
  annotation_col = mb_srat[[]]["subtype"]
)

# ------------------------------------------------------------------------------
# SingleR with combined RL lineage and UBC subcluster labels as reference

message("***Projecting MB tumour cells on RL lineage and UBC subclusters***")

subdir <- file.path(arg_list$out_dir, "merged_RL_with_UBC")
if (!dir.exists(subdir)) {
  dir.create(subdir, recursive = TRUE)
}

# combine RL lineage labels with UBC subcluster labels to create the reference
# annotations
rl_lineage <- rl_srat[[]][, arg_list$rl_srat[2]] %>%
  as.character() %>%
  setNames(nm = rownames(rl_srat[[]]))
ubc_clust <- ubc_srat[[]][, arg_list$ubc_srat[2]] %>% 
  paste0("UBC_", .) %>%
  setNames(nm = rownames(ubc_srat[[]]))
merged_labels <- replace(
  x = rl_lineage,
  list = names(ubc_clust),
  values = ubc_clust
)
merged_labels <- merged_labels[colnames(rl_mat)] %>% unname()

# run SingleR
annot_merged <- SingleR(
  test = mb_mat,
  ref = rl_mat,
  labels = merged_labels,
  de.method = "wilcox",
  num.threads = arg_list$num_threads
)
message(sprintf(
  "***%s labels were pruned***",
  sum(is.na(annot_merged$pruned.labels))
))

saveRDS(annot_merged, file = file.path(subdir, "prediction_results.rds"))

# add labels back to metadata
mb_srat$merged_rl_ubc_pred <- annot_merged[rownames(mb_srat[[]]), "pruned.labels"]

# save plots
make_plots(
  srat = mb_srat,
  srat_pred_col = "merged_rl_ubc_pred",
  path = subdir,
  prediction_results = annot_merged,
  annotation_col = mb_srat[[]]["subtype"]
)

# ------------------------------------------------------------------------------
# save tumour metadata with predicted cell types

saveRDS(
  object = mb_srat[[]],
  file = file.path(
    arg_list$out_dir,
    paste0(tools::file_path_sans_ext(basename(arg_list$mb_srat[1])), "_metadata.rds")
  )
)

# ------------------------------------------------------------------------------

message("\n***SESSION INFO***\n")
print(sessionInfo())

