# ==============================================================================
# Make a Venn diagram showing the human-specific cells that are shared across
# integration methods.
# ==============================================================================

# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(argparse)
library(tidyverse)
library(Seurat)
library(ggVennDiagram)

# parse command line arguments
parser <- ArgumentParser()
parser$add_argument("--date_dir", default = NULL)

args <- parser$parse_args()

# set data and output directories
root_dir <- "CBL_scRNAseq"
out_dir <- file.path("", root_dir, "results/integrated")

if (is.null(args$date_dir)) {
  date_dir <- file.path(out_dir, format(Sys.Date(), "%Y%m%d"))
} else {
  date_dir <- args$date_dir
}

if (!dir.exists(date_dir)) {
  dir.create(date_dir)
}

message(sprintf("Saving results to: %s", date_dir))

# import functions
source(file.path("", root_dir, "software/utilities/cell_labelling.R"))

# ------------------------------------------------------------------------------
# load Seurat objects

cca_srat <- readRDS(file.path(out_dir, "vladoiu_liam_RL.rds"))
rpca_srat <- readRDS(file.path(out_dir, "20221002/kanchors80/vladoiu_liam_RL_rpca.rds"))
harmony_srat <- readRDS(file.path(out_dir, "20221003/vladoiu_liam_RL_harmony.rds"))

# ------------------------------------------------------------------------------
# identify human-specific UBCs and make Venn diagram

human_spec_ubcs <- list(
  CCA = WhichCells(cca_srat,
                   # cluster 19 and 20 are human-specific UBCs
                   expression = seurat_clusters %in% c(19, 20) &
                     # remove stray mouse cells from cluster
                     species == "human"), 
  RPCA = WhichCells(rpca_srat,
                    # clusters 19 and 22 are human-specific UBCs
                    expression = seurat_clusters %in% c(19, 22) &
                      # remove stray mouse cells from cluster
                      species == "human"),
  harmony = WhichCells(harmony_srat,
                       # cluster 16 is human-specific UBCs
                       expression = seurat_clusters %in% c(16) &
                         # remove stray mouse cells from cluster
                         species == "human")
)

plt <- ggVennDiagram(human_spec_ubcs, label_alpha = 0, edge_size = 0) + 
  scale_x_continuous(expand = expansion(mult = 0.2)) + 
  scale_fill_distiller(palette = "Blues", direction = 1)
ggsave("human_specific_ubcs_venn.png", plot = plt, path = date_dir, width = 6, height = 4, units = "in", dpi = 600)
