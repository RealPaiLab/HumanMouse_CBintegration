# ==============================================================================
# Run GSEA with g:Profiler on non-orthologous genes to see what data is being
# lost during integration.
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
  # CSV of human/mouse orthologs
  "--orth_csv",
  default = NULL,
  required = TRUE
)
parser$add_argument(
  # custom GMT gene set
  "--organism",
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
    "--orth_csv", "/.mounts/labs/pailab/private/icheong/CBL_scRNAseq/results/integrated/hgnc_mgi_orth_genes_20240402.csv",
    "--out_dir", "/.mounts/labs/pailab/private/icheong/CBL_scRNAseq/results/integrated/20241005"
  ))
} else {
  arg_list <- parser$parse_args()
}

message(sprintf("Saving files to %s", getwd()))
if (!dir.exists(arg_list$out_dir)) {
  dir.create(arg_list$out_dir, recursive = TRUE)
}

# load functions
source("/.mounts/labs/pailab/private/icheong/CBL_scRNAseq/software/utilities/gprofiler2_helpers.R")

# ------------------------------------------------------------------------------
# load human/mouse orthologous genes from file

message(sprintf("***Reading in human/mouse orthologs from %s***", arg_list$orth_csv))
orth_genes <- read_csv(arg_list$orth_csv)

# remove any orthologs that aren't one-to-one; taken from Leo's code
# (/.mounts/labs/pailab/private/icheong/CBL_scRNAseq/scrnaseq_Leo/utilities/mouse_to_human_genes.R)
orth_genes <- orth_genes[
  !(duplicated(orth_genes$MGI.symbol) | duplicated(orth_genes$MGI.symbol, fromLast = TRUE))
  & !(duplicated(orth_genes$HGNC.symbol) | duplicated(orth_genes$HGNC.symbol, fromLast = TRUE))
  & orth_genes$MGI.symbol != "Pisd",
]

# save one-to-one orthologs
write_csv(
  x = orth_genes,
  file = file.path(arg_list$out_dir, "human_mouse_one_to_one_orthologs.csv")
)

# ------------------------------------------------------------------------------
# path to individual datasets

srat_qs <- read_csv(
  "/.mounts/labs/pailab/private/icheong/CBL_scRNAseq/scrnaseq_Leo/integrations/dataset_list.csv"
) %>%
  dplyr::select(-dataset_location_rds) %>%
  dplyr::filter(
    str_detect(dataset_name, "_full_cerebellum_"),
    !str_detect(dataset_name, "Luo")
  ) %>%
  mutate(
    dataset_name = tolower(str_remove(dataset_name, "full_cerebellum_"))
  ) %>%
  tibble::deframe()

# ------------------------------------------------------------------------------
# run gost on human genes that don't have one-to-one orthologs

# load human datasets
srat <- map(
  .x = c("aldinger_human", "sepp_human"),
  .f = \(dataset) {
    message(sprintf("***Reading in %s from %s***", dataset, srat_qs[dataset]))
    qs::qread(srat_qs[dataset])
  }
) %>%
  setNames(nm = c("aldinger_human", "sepp_human"))

# get human genes
hsap_gene <- union(rownames(srat$aldinger_human), rownames(srat$sepp_human))
write_lines(
  x = hsap_gene,
  file = file.path(arg_list$out_dir, "aldinger_sepp_human_genes.txt")
)

# get non-orthologous human genes
nonorth_human <- setdiff(hsap_gene, orth_genes$HGNC.symbol)
message(sprintf("***There are %s non-orthologous human genes***", length(nonorth_human)))

gost_res <- gost(
  query = nonorth_human,
  organism = arg_list$organism,
  significant = TRUE,
  evcodes = TRUE,
  correction_method = "fdr",
  custom_bg = hsap_gene
)

write_csv(
  x = gost_res$result,
  file = file.path(arg_list$out_dir, "nonorth_human.csv")
)
write_rds(
  x = gost_res,
  file = file.path(arg_list$out_dir, "nonorth_human.rds")
)

# save an EnrichmentMap format
gost_res2gem(
  gost_res = gost_res,
  phenotype = "+1"
) %>% 
  write_gem(file = file.path(arg_list$out_dir, "nonorth_human.gem.txt"))

# save plot
plt <- gostplot(gost_res)
htmlwidgets::saveWidget(
  widget = plt,
  file = file.path(arg_list$out_dir, "nonorth_human-manhattan.html"),
  title = "non-orthologous human genes"
)

# ------------------------------------------------------------------------------
# run gost on mouse genes that don't have one-to-one orthologs

#TODO same as above but for mouse

# ------------------------------------------------------------------------------

message("\n***SESSION INFO***\n")
print(sessionInfo())

