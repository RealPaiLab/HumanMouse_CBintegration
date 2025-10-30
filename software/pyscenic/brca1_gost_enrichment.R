# ==============================================================================
# Run GSEA with g:Profiler on BRCA1 regulons to see what pathways are active.
# ==============================================================================

# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(argparse)
library(tidyverse)
library(Seurat)
library(biomaRt)
library(gprofiler2)

# parse command line arguments
parser <- ArgumentParser()

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
    "--out_dir", "/CBL_scRNAseq/results/integrated/20231006"
  ))
} else {
  arg_list <- parser$parse_args()
}

message(sprintf("Saving files to %s", arg_list$out_dir))
if (!dir.exists(arg_list$out_dir)) {
  dir.create(arg_list$out_dir, recursive = TRUE)
}

# load functions
source("/.mounts/labs/pailab/private/icheong/CBL_scRNAseq/software/utilities/gprofiler2_helpers.R")

# ------------------------------------------------------------------------------
# load regulons and background genes

# BRCA1 in normal development
rl_brca1 <- read_csv("/.mounts/labs/pailab/private/icheong/CBL_scRNAseq/results/pyscenic/20241017/BRCA1.csv")
# background genes
rl_srat <- read_rds("/.mounts/labs/pailab/private/icheong/CBL_scRNAseq-archived/data/human/Aldinger/glutamatergic_dev_Liam.RDS")
rl_genes <- rownames(rl_srat)

# BRCA1 in tumours
mb_brca1 <- read_csv("/.mounts/labs/pailab/private/icheong/CBL_scRNAseq/results/pyscenic/20231019/BRCA1.csv")
# background genes
mb_srat <- read_rds("/.mounts/labs/pailab/private/icheong/CBL_scRNAseq/results/tumour/Vladoiu/20230510/mb_mnn.rds")
mb_genes <- rownames(mb_srat)

# ------------------------------------------------------------------------------
# run gost

regulon <- list(rl_brca1, mb_brca1)
bg <- list(rl_genes, mb_genes)
subdir <- list("normal", "tumour")

pwalk(
  .l = list(regulon, bg, subdir),
  .f = \(
    regulon,
    bg,
    subdir,
    organism = arg_list$organism,
    out_dir = arg_list$out_dir
  ) {
    out_dir <- file.path(out_dir, subdir)
    if (!dir.exists(out_dir)) {
      dir.create(out_dir, recursive = TRUE)
    }
    
    # run gost
    message(sprintf("\n***Running g:GOSt on %s***", subdir))
    gost_res <- gost(
      query = regulon$gene,
      organism = organism,
      significant = TRUE,
      evcodes = TRUE,
      correction_method = "fdr",
      custom_bg = bg
    )
    
    # save output
    write_csv(
      x = gost_res$result,
      file = file.path(out_dir, paste0(subdir, ".gost.csv"))
    )
    write_rds(
      x = gost_res,
      file = file.path(out_dir, paste0(subdir, ".gost.rds"))
    )
    
    # save as EnrichmentMap format
    gost_gem <- gost_res2gem(gost_res = gost_res, phenotype = "+1") %>% 
      write_gem(file.path(out_dir, paste0(subdir, ".gem.txt")))
    
    # save enriched pathways as `.gem` file for Cytoscape
    gost_res2gem(
      gost_res = gost_res,
      # min_max_terms = c(10, 250),
      phenotype = "+1"
    ) %>% 
      write_gem(file.path(out_dir, paste0(subdir, "_min10_max250.gem.txt")))
    
    # save plot
    plt <- gostplot(gost_res)
    htmlwidgets::saveWidget(
      widget = plt,
      file = file.path(out_dir, paste0(subdir, "-manhattan.html")),
      title = paste(subdir, "BRCA1 regulon")
    )
  }
)

# ------------------------------------------------------------------------------

message("\n***SESSION INFO***\n")
print(sessionInfo())

