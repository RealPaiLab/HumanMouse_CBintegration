# ==============================================================================
# Run GSEA with g:Profiler on edgeR results.
# ==============================================================================

# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(argparse)
library(tidyverse)
library(gprofiler2)

# parse command line arguments
parser <- ArgumentParser()

parser$add_argument(
  # input CSV files derived from a `TopTags` object
  "--genes_csv",
  action = "extend",
  nargs = "+", # need at least one file
  default = NULL,
  required = TRUE
)
parser$add_argument(
  # organism to use in `gost()`
  "--organism",
  default = "gp__CEZi_x2oO_5Lg",
  help = "passed to the `organism` parameter in `gost()`"
)
parser$add_argument(
  # output directory
  "--out_dir",
  default = NULL,
  required = TRUE
)

if (interactive()) {
  arg_list <- parser$parse_args(c(
    "--genes_csv", dir("/CBL_scRNAseq/results/species_dge/20230906/",
                       pattern = "*.csv", full.names = TRUE),
    "--out_dir", "/CBL_scRNAseq/results/species_dge/20230914"
  ))
} else {
  arg_list <- parser$parse_args()
}

message(sprintf("Saving files to %s", getwd()))
if (!dir.exists(arg_list$out_dir)) {
  dir.create(arg_list$out_dir, recursive = TRUE)
}

# load functions
source("/CBL_scRNAseq/software/utilities/gprofiler2_helpers.R")

# ------------------------------------------------------------------------------
# load data

de_genes <- map(
  .x = arg_list$genes_csv,
  .f = read_csv
)
names(de_genes) <- tools::file_path_sans_ext(basename(arg_list$genes_csv))

# ------------------------------------------------------------------------------
# run gost for each set of differentially expressed genes

walk2(
  .x = de_genes,
  .y = names(de_genes),
  .f = \(x, y) {
    # extract the up- and down-regulated genes
    down_genes <- x$gene[x$logFC < 0 & x$PValue < 0.05]
    up_genes <- x$gene[x$logFC > 0 & x$PValue < 0.05]
    
    query <- list(down_genes, up_genes)
    names(query) <- paste0(y, c(" - down", " - up"))
    
    # run gost
    message(sprintf("\n***Running g:GOSt on %s***", y))
    gost_res <- gost(
      query = query,
      organism = arg_list$organism,
      significant = TRUE, # return only significant results
      evcodes = TRUE, # required to get intersection of genes
      correction_method = "fdr",
      custom_bg = x$gene
    )
    
    # save results
    write_csv(
      x = gost_res$result,
      file = file.path(arg_list$out_dir, paste0(y, ".gost.csv"))
    )
    write_rds(
      x = gost_res,
      file = file.path(arg_list$out_dir, paste0(y, ".gost.rds"))
    )
    
    # save as EnrichmentMap format
    gost_res2gem(
      gost_res = gost_res,
      phenotype = list(down = names(query)[1], up = names(query)[2])
    ) %>% 
      write_gem(file = file.path(arg_list$out_dir, paste0(y, ".gem.txt")))
    
    # make Manhattan plot and save
    plt <- gostplot(gost_res)
    htmlwidgets::saveWidget(
      widget = plt,
      file = file.path(arg_list$out_dir, paste0(y, "-manhattan.html")),
      title = y
    )
  }
)

# ------------------------------------------------------------------------------

message("\n***SESSION INFO***\n")
print(sessionInfo())
