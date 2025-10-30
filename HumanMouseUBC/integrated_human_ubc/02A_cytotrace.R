# ==============================================================================
# Run CytoTRACE on integrated human UBC subclusters.
# ==============================================================================

library(argparse)
library(tidyverse)
library(Seurat)

# set python for CytoTRACE
reticulate::use_condaenv("monocle3")
library(CytoTRACE)

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
  # number of cores to use
  "--num_threads",
  default = 1,
  required = FALSE
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
    "--out_dir", "/.mounts/labs/pailab/private/projects/HumanMouseUBC/integrated_human_ubc/20250312/"
  ))
} else {
  arg_list <- parser$parse_args()
}

message(sprintf("Saving files to %s", arg_list$out_dir))
if (!dir.exists(arg_list$out_dir)) {
  dir.create(arg_list$out_dir, recursive = TRUE)
}

# ------------------------------------------------------------------------------
# functions

#' Wrapper function to run CytoTRACE on a Seurat object.
#' 
#' @param srat Seurat object.
#' @param num_threads Number of threads/cores to use for running CytoTRACE
#'   (passed to `ncores`).
#' @param batch Metadata column name to run CytoTRACE on heterogenous
#'   batches/datasets. If `NULL`, run the default `CytoTRACE()` function without
#'   correcting for any batch effects. Otherwise, `iCytotrace()` is run which
#'   uses `scanorama` to normalize and merge heterogenous datasets; see
#'   CytoTRACE documentation for more details.
#' @param out_dir Path to save the CytoTRACE output to.
#' 
#' @return Output of running `CytoTRACE`/`iCytoTRACE`.
#' 
run_ct <- function(
  srat,
  num_threads,
  batch = NULL,
  out_dir
) {
  if (is.null(batch)) {
    # don't split Seurat object (run all cells together)

    # get raw counts (CytoTRACE requires dense matrix)
    raw_counts <- GetAssayData(srat, assay = "RNA", slot = "counts") %>%
      as.matrix()

    message("***Running `CytoTRACE`***")
    t0 <- Sys.time()
    ct_res <- CytoTRACE(mat = raw_counts, ncores = num_threads)
    t1 <- Sys.time()
    message(sprintf(
      "***CytoTRACE runtime: %.02f min***",
      difftime(t1, t0, units = "mins")
    ))
  } else {
    # split Seurat object by `batch`

    # get raw counts split by batch
    srat_batched <- SplitObject(srat, split.by = batch)
    raw_counts <- map(
      .x = srat_batched,
      .f = \(x) {
        GetAssayData(x, assay = "RNA", slot = "counts") %>%
          as.matrix()
      }
    ) %>%
      unname()

    message("***Running `iCytoTRACE`***")
    t0 <- Sys.time()
    ct_res <- iCytoTRACE(datasets = raw_counts, ncores = num_threads)
    t1 <- Sys.time()
    message(sprintf(
      "***iCytoTRACE runtime: %.02f min***",
      difftime(t1, t0, units = "mins")
    ))
  }

  # save CytoTRACE results
  saveRDS(ct_res, file.path(out_dir, "cytotrace_output.rds"))

  return(ct_res)
}

#' Make violin plots of CytoTRACE scores coloured/split by different groups.
#' 
#' @param srat Seurat object.
#' @param groups At least one of `c("age", "dataset", "cluster")` and/or another
#'   column present in the Seurat metadata.
#' @param out_dir Path to save the plots to.
#' 
make_violin <- function(
  srat,
  groups,
  out_dir
) {
  walk(
    .x = groups,
    .f = \(group_by) {
      # get Seurat metadata
      md <- srat[[]]

      # filename for plot
      filename <- sprintf("vln_%s.png", group_by)

      # make violin plot
      group_by <- switch(
        group_by,
        dataset = "dataset_name",
        cluster = arg_list$cluster,
        group_by
      )
      .plt <- ggplot(md, aes(x = .data[[group_by]], y = cytotrace)) + 
        geom_violin(aes(fill = .data[[group_by]])) + 
        geom_jitter(width = 0.25, size = 0.1) + 
        labs(x = group_by, y = "CytoTRACE score\n(higher = more stem-like)") + 
        theme_classic() + 
        theme(
          axis.text = element_text(colour = "black"),
          axis.ticks = element_line(colour = "black"),
          legend.position = "none"
        )

      # save plot
      ggsave(
        plot = .plt,
        filename = filename,
        path = out_dir,
        width = length(unique(md[[group_by]])) + 1,
        height = 6,
        units = "in",
        dpi = 600
      )
    }
  )
}

# ------------------------------------------------------------------------------
# load Seurat object

srat <- readRDS(arg_list$srat_file)
srat$dataset_name <- str_remove(srat$dataset_name, "_full_cerebellum_human")
srat$age <- factor(srat$age) %>%
  fct_relevel(\(x) {str_sort(x, numeric = TRUE)})

# ------------------------------------------------------------------------------
# run CytoTRACE and make plots

pwalk(
  .l = list(
    batch = list(NULL, "dataset_name"),
    subdir = list("without_correction", "with_dataset_correction")
  ),
  .f = \(batch, subdir) {
    # set subdirectory for results
    subdir <- file.path(arg_list$out_dir, subdir)
    if (!dir.exists(subdir)) {
      dir.create(subdir, recursive = TRUE)
    }

    # run CytoTRACE
    ct_res <- run_ct(
      srat = srat,
      num_threads = arg_list$num_threads,
      batch = batch,
      out_dir = subdir
    )

    # add CytoTRACE scores to Seurat object metadata
    srat$cytotrace <- ct_res$CytoTRACE[rownames(srat[[]])]    

    # get order of clusters based on CytoTRACE (from least to most differentiated)
    ubc_clust_order <- srat[[]] %>%
      select(all_of(arg_list$cluster), cytotrace) %>%
      group_by(.data[[arg_list$cluster]]) %>%
      summarise(median = median(cytotrace)) %>%
      arrange(desc(median)) %>%
      pull(.data[[arg_list$cluster]])

    # order cluster levels
    srat@meta.data[[arg_list$cluster]] <- factor(
      srat[[]][[arg_list$cluster]],
      levels = ubc_clust_order
    )

    # >>> plotting >>>

    # UMAP showing CytoTRACE scores
    .plt <- FeaturePlot(
      srat,
      features = "cytotrace",
      reduction = "umap"
    ) + 
      labs(title = "CytoTRACE scores") + 
      scale_colour_viridis_c(direction = -1)
    ggsave(
      plot = .plt,
      filename = "umap_cytotrace.png",
      path = subdir,
      width = 5.5,
      height = 5,
      units = "in",
      dpi = 600
    )

    # violin plots of cytotrace results by age/dataset/cluster
    make_violin(
      srat = srat,
      groups = c("age", "dataset", "cluster"),
      out_dir = subdir
    )

    # <<<
  }
)

# ------------------------------------------------------------------------------
# additional UMAP plots

# UMAPs of CytoTRACE results (show cluster as well)
walk(
  .x = c("cluster", "dataset"),
  .f = \(group_by) {
    # filename for plot
    filename <- sprintf("umap_%s.png", group_by)

    # make UMAP
    group_by <- switch(
      group_by,
      cluster = arg_list$cluster,
      dataset = "dataset_name",
    )
    .plt <- DimPlot(
      srat,
      reduction = "umap",
      group.by = group_by,
      label = TRUE,
      repel = TRUE
    ) + 
      NoLegend()

    # save plot
    ggsave(
      plot = .plt,
      filename = filename,
      path = arg_list$out_dir,
      width = 5,
      height = 5,
      units = "in",
      dpi = 600
    )
  }
)

# ------------------------------------------------------------------------------

message("\n***SESSION INFO***\n")
print(sessionInfo())

