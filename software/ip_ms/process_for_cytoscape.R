# ==============================================================================
# Add known molecular interaction data and known DNA-binding proteins to BRCA1
# IP-MS results for Cytoscape visualization.
# ==============================================================================

# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(argparse)
library(tidyverse)

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
    "--out_dir", "/CBL_scRNAseq/results/ip_ms/20240109/"
  ))
} else {
  arg_list <- parser$parse_args()
}

# path to input Excel file
filename <- "/isilon/CBL_scRNAseq-archived/data/ip_ms/BRCA1_Kaiso_interaction_results_20240102.xlsx"
# name of sheet from input Excel file
sheetname <- "BRCA1_SAINTID6835"
# path to CORUM database file
corum_file <- "/isilon/CBL_scRNAseq-archived/data/src/molecular-association/CORUM/CORUM-download-2022_09_12.xlsx"
# path to human transcription factors database (with DNA-binding domains)
human_tfs_file <- "/isilon/CBL_scRNAseq-archived/data/src/transcription-factors/human-tfs-lambert/full_database_v1.01.csv"
# output directory
out_dir <- arg_list$out_dir

# load functions
source("/CBL_scRNAseq/software/ip_ms/load_ip_ms_data.R")

# ------------------------------------------------------------------------------
# find proteins that form complexes with BRCA1

# CORUM database + filtering
corum_db <- readxl::read_excel(path = corum_file) %>% 
  # filter for human complexes containing BRCA1
  filter(
    Organism == "Human",
    grepl("P38398", `subunits(UniProt IDs)`) # BRCA1 UniProt ID
  )

# get list of proteins
complex_proteins <- corum_db$`subunits(UniProt IDs)` %>% 
  str_split(pattern = ";") %>% 
  unlist() %>% 
  unique()

write_lines(
  x = complex_proteins,
  file = file.path(out_dir, "proteins_complexing_with_brca1.txt")
)

# ------------------------------------------------------------------------------
# find DNA-binding proteins in Lambert's database of human TFs

# load databse + filtering
dna_binding <- read_csv(file = human_tfs_file) %>% 
  filter(`Binding mode` != "Not a DNA binding protein") %>% 
  select("HGNC symbol") %>% 
  unlist() %>% 
  unique()

write_lines(
  x = dna_binding,
  file = file.path(out_dir, "dna_binding_proteins.txt")
)

# ------------------------------------------------------------------------------
# add column to IP-MS data and save as new file

dat <- load_ip_ms_data(filename, sheetname) %>% 
  mutate(
    `complexes with BRCA1 (CORUM)` = if_else(
      Prey %in% complex_proteins,
      "yes",
      "no"
    ),
    `binds DNA` = if_else(
      PreyGene %in% dna_binding,
      "yes",
      "no"
    )
  )

write_csv(
  x = dat,
  file = file.path(out_dir, "brca1_ip_ms.csv")
)

# ------------------------------------------------------------------------------

message("\n***SESSION INFO***\n")
print(sessionInfo())

