# ==============================================================================
# Retrieve Ensembl gene IDs and HGNC symbols from biomart for the Sepp dataset.
# ==============================================================================

# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(argparse)
library(SingleCellExperiment)
library(biomaRt)

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
    "--out_dir", "/CBL_scRNAseq/results/human/Sepp/20240329/"
  ))
} else {
  arg_list <- parser$parse_args()
}

message(sprintf("***Saving files to %s***", arg_list$out_dir))
if (!dir.exists(arg_list$out_dir)) {
  dir.create(arg_list$out_dir, recursive = TRUE)
}

# ------------------------------------------------------------------------------
# load SCE object and get genes

sce <- readRDS("/isilon/CBL_scRNAseq-archived/data/src/neurodev-genomics/scRNAseq/Sepp_2024/hum_sce_final.rds")

ensg_genes <- rownames(sce)

# ------------------------------------------------------------------------------
# get HGNC genes

# load biomart
bm <- useEnsembl(
  biomart = "genes",
  dataset = "hsapiens_gene_ensembl",
  host = "https://jan2024.archive.ensembl.org",
  mirror = "useast"
)

# query biomart
gene_table <- getBM(
  attributes = c("ensembl_gene_id", "hgnc_symbol"),
  filters = "ensembl_gene_id",
  values = ensg_genes,
  mart = bm,
  uniqueRows = FALSE
)

# ------------------------------------------------------------------------------
# save output

readr::write_csv(
  x = gene_table,
  file = file.path(arg_list$out_dir, "ensg_hgnc_ids.csv")
)

# ------------------------------------------------------------------------------

message("\n***SESSION INFO***\n")
print(sessionInfo())
