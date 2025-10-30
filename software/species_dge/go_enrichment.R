# ==============================================================================
# GO enrichment of species-specific genes.
# ==============================================================================

# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(argparse)
library(tidyverse)
library(edgeR)
library(AnnotationDbi)
library(org.Hs.eg.db)

# parse command line arguments
parser <- ArgumentParser()

parser$add_argument(
  # input CSV files
  "--genes_csv",
  action = "extend",
  nargs = "+", # need at least one file
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
    "--genes_csv", dir("/CBL_scRNAseq/results/species_dge/20230818/",
                       pattern = "*.csv", full.names = TRUE),
    "--out_dir", "/CBL_scRNAseq/results/species_dge/20230824/"
  ))
} else {
  arg_list <- parser$parse_args()
}

message(sprintf("Saving files to %s", getwd()))
if (!dir.exists(arg_list$out_dir)) {
  dir.create(arg_list$out_dir, recursive = TRUE)
}

# ------------------------------------------------------------------------------
# load data

cell_types <- tools::file_path_sans_ext(basename(arg_list$genes_csv))
all_de_genes <- map(
  .x = arg_list$genes_csv,
  .f = read_csv
) %>% 
  set_names(cell_types)

# ------------------------------------------------------------------------------
# GO enrichment

#' Find enriched GO terms in a gene list
#'
#' @param sig_genes Differentially expressed genes.
#' @param uni_genes Universe of genes for background.
#'
#' @return Enriched GO terms.
#'
find_enr_go <- function(
  sig_genes,
  ctrl_genes
) {
  # convert to Entrez Gene IDs
  sig_genes <- mapIds(
    x = org.Hs.eg.db,
    keys = sig_genes,
    column = "ENTREZID",
    keytype = "SYMBOL"
  )
  ctrl_genes <- mapIds(
    x = org.Hs.eg.db,
    keys = ctrl_genes,
    column = "ENTREZID",
    keytype = "SYMBOL"
  )
  
  # GO enrichment
  enr_go <- goana(
    de = sig_genes,
    universe = ctrl_genes
  ) %>% 
    topGO(number = Inf)
  
  return(enr_go)
}

# for human-specific genes (logFC < -0.3)
hs_sig_genes <- map(
  .x = all_de_genes,
  .f = \(x) {
    sig_genes <- x[x$logFC < -0.3 & x$FDR < 0.05, ]$gene
    return(sig_genes)
  }
)
hs_ctrl_genes <- map(
  .x = all_de_genes,
  .f = \(x) {x$gene}
)

# human-specific GO for each cell type
message("Running GO term enrichment for each cell type")
hs_enr_go <- map2(
  .x = hs_sig_genes,
  .y = hs_ctrl_genes,
  .f = \(x, y) {
    enr_go <- find_enr_go(x, y) %>% 
      rownames_to_column(var = "GO.ID")
    return(enr_go)
  }
)

# save to file
walk2(
  .x = hs_enr_go,
  .y = names(hs_enr_go),
  .f = \(x, y) {
    write_csv(
      x = x, 
      file = file.path(arg_list$out_dir, paste0(y, ".hs_go.csv")))
  }
)

# ------------------------------------------------------------------------------

message("\n***SESSION INFO***\n")
print(sessionInfo())

