# ==============================================================================
# Retrieve human and mouse orthologous genes from BioMart. Run this script to
# get all human-mouse orthologs from BioMart and save it as a CSV file (date
# included in the filename).
# ==============================================================================

# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(biomaRt)

# set output directory
out_dir <- "/CBL_scRNAseq/results/integrated/"

message(sprintf("***Saving files to %s***", out_dir))
if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

# load functions
source("/CBL_scRNAseq/software/utilities/convert_genes.R")

# ------------------------------------------------------------------------------
# load biomart

human_bm <- load_human_bm()
mouse_bm <- load_mouse_bm()

# ------------------------------------------------------------------------------
# get orthologs

# mouse orthologs of human genes
message("Retrieving mouse orthologs of human genes")
mouse_genes <- convert_h2m(
  human_bm = human_bm,
  mouse_bm = mouse_bm
)

# human orthologs of mouse genes
message("Retrieving human orthologs of mouse genes")
human_genes <- convert_m2h(
  mouse_bm = mouse_bm,
  human_bm = human_bm
)

# combine both tables (rbind automatically matches by column name)
all_orth_genes <- rbind(mouse_genes, human_genes) |>
  dplyr::distinct() |> 
  dplyr::filter(HGNC.symbol != "" & MGI.symbol != "")

# save output
readr::write_csv(
  x = all_orth_genes,
  file = file.path(out_dir, paste0("hgnc_mgi_orth_genes_", format(Sys.time(), "%Y%m%d"), ".csv"))
)
