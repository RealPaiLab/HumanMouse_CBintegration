# ==============================================================================
# Import ICGC MB mutation data and convert into a suitable format for Oncoprint
# ==============================================================================

# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(argparse)
library(tidyverse)
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
    "--out_dir", "/CBL_scRNAseq/results/tumour/ICGC/20231102"
  ))
} else {
  arg_list <- parser$parse_args()
}

# load functions
# source("/CBL_scRNAseq/software/utilities/convert_genes.R")

# ------------------------------------------------------------------------------
# load data

peme <- read_tsv("/isilon/CBL_scRNAseq-archived/data/src/ICGC/PEME-CA/simple_somatic_mutation.open.PEME-CA.tsv.gz") %>% 
  dplyr::select(c(
    icgc_donor_id,
    # project_code,
    # submitted_sample_id,
    mutation_type,
    gene_affected
  )) %>% 
  unique() %>% 
  filter(
    !is.na(gene_affected)
  )

# load biomaRt db; need GRCh37 cuz that's what was used in the PEME-CA data
message("***Loading BioMart***")
human_bm <- useEnsembl(
  biomart = "genes",
  dataset = "hsapiens_gene_ensembl",
  GRCh = 37,
  verbose = TRUE
)
gene_conversion <- getBM(
  attributes = c("ensembl_gene_id", "hgnc_symbol"),
  mart = human_bm,
  filters = "ensembl_gene_id",
  values = peme$gene_affected
)
write_csv(gene_conversion, file = file.path(arg_list$out_dir, "ensembl_to_hgnc.csv"))

# ------------------------------------------------------------------------------
# convert to oncoprint matrix format

# pivot wider to matrix format
peme_mat <- peme %>% 
  mutate(
    mutation_type = case_when(
      mutation_type == "single base substitution" ~ "snv",
      mutation_type == "deletion of <=200bp" ~ "indel <=200",
      mutation_type == "insertion of <=200bp" ~ "indel <=200"
    )
  ) %>% 
  pivot_wider(
    names_from = icgc_donor_id,
    values_from = mutation_type,
    values_fn = \(x) {
      # get only unique values and collapse
      x <- unique(x) %>% 
        sort() %>% 
        paste0(collapse = ";")
      return(x)
    },
    values_fill = ""
  ) %>% 
  left_join(
    # convert ensembl ID to gene symbol
    y = gene_conversion,
    by = c("gene_affected" = "ensembl_gene_id")
  ) %>% 
  filter(hgnc_symbol != "" & !duplicated(hgnc_symbol)) %>% 
  dplyr::select(-gene_affected) %>% 
  column_to_rownames("hgnc_symbol") %>% 
  as.matrix()

# save as CSV
write.csv(
  peme_mat,
  file = file.path(arg_list$out_dir, "simple_somatic_mutation.PEME-CA.csv")
)


#TODO filter data for oncoprint; take only genes with highest frequency


# ------------------------------------------------------------------------------

message("\n***SESSION INFO***\n")
print(sessionInfo())

