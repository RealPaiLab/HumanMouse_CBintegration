# ==============================================================================
# Check diferentially expressed genes in UBC subclusters for enriched pathways.
# ==============================================================================

library(argparse)
library(tidyverse)
library(gprofiler2)

# parse command line arguments
parser <- ArgumentParser()

parser$add_argument(
  # differential gene expression results
  "--de_genes",
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
    "--de_genes", "/.mounts/labs/pailab/private/projects/HumanMouseUBC/integrated_human_ubc/20250331/regress_sex_dataset/de_genes.tsv",
    "--out_dir", "/.mounts/labs/pailab/private/projects/HumanMouseUBC/integrated_human_ubc/PathwayEnrichment"
  ))
} else {
  arg_list <- parser$parse_args()
}

dt <- format(Sys.Date(),"%y%m%d")
arg_list$out_dir <- file.path(arg_list$out_dir, dt)

message(sprintf("Saving files to %s", getwd()))
if (!dir.exists(arg_list$out_dir)) {
  dir.create(arg_list$out_dir, recursive = TRUE)
}

# ------------------------------------------------------------------------------
# functions

source("/.mounts/labs/pailab/private/icheong/CBL_scRNAseq/software/utilities/gprofiler2_helpers.R")

get_upreg_genes <- function(
  de_res
) {
  upreg_genes <- de_res %>%
    filter(
      if_all(ends_with("avg_log2FC"), ~ .x > 0),
      if_all(ends_with("p_val") & !ends_with("minimump_p_val"), ~ .x < 0.05)
    ) %>%
    pull(gene)
  return(upreg_genes)
}

# ------------------------------------------------------------------------------
# load differential gene expression results

message(sprintf(
  "***Reading in differential gene expression results from %s***",
  arg_list$de_genes
))
de_genes <- read_tsv(arg_list$de_genes)

# ------------------------------------------------------------------------------
# run pathway enrichment analysis with g:Profiler

# get UBC subclusters
ubc_subcluster <- unique(de_genes$ubc_subcluster)

# run g:Profiler for upregulated genes in each UBC subcluster
walk(
  .x = ubc_subcluster,
  .f = \(clust) {
    # create cluster subdirectory
    subdir <- file.path(arg_list$out_dir, clust)
    if (!dir.exists(subdir)) {
      dir.create(subdir, recursive = TRUE)
    }

    # get results for the cluster
    clust_res <- de_genes[de_genes$ubc_subcluster == clust, ]
    upreg_genes <- get_upreg_genes(de_res = clust_res)
    bg_genes <- clust_res$gene

    message(sprintf(
      "***%s/%s genes are upregulated in both datasets (p < 0.05)***",
      length(upreg_genes),
      length(bg_genes)
    ))

    # save upregulated and background genes
    write_lines(upreg_genes, file.path(subdir, "upreg_genes.txt"))
    write_lines(bg_genes, file.path(subdir, "bg_genes.txt"))

    cat("about to run g:Profiler with custom pathway file\n")

    # run pathway enrichment analysis (from `gprofiler2_helpers.R`)
    gost_res <- run_gost(
      query = upreg_genes,
      organism = "gp__OUFl_gpIE_RyE",
      significant = TRUE,
      evcodes = TRUE,
      correction_method = "fdr",
      custom_bg = bg_genes,
      filename = file.path(subdir, "enriched_pathways")
    )

    # save results in Generic Enrichment Map format for downstream Cytoscape
    # (functions from `gprofiler2_helpers.R`)
    if (!is.null(gost_res)) {
      gost_res2gem(gost_res = gost_res, phenotype = "+1") %>%
        write_gem(file = file.path(subdir, "enriched_pathways.gem.txt"))
    }
  }
)

# ------------------------------------------------------------------------------

message("\n***SESSION INFO***\n")
print(sessionInfo())