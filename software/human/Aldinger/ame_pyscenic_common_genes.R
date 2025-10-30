# ==============================================================================
# Get intersection of TFs from MEMES and pySCENIC analysis.
# ==============================================================================

# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(argparse)
library(tidyverse)
library(ggVennDiagram)
library(patchwork)

# parse command line arguments
parser <- ArgumentParser()

parser$add_argument(
  # path to CSV of enriched motifs from AME/MEMES
  "--ame_csv",
  default = NULL,
  required = TRUE
)
parser$add_argument(
  # path to CSV of RSS scores
  "--rss_csv",
  default = NULL,
  required = TRUE
)
parser$add_argument(
  # cell type cluster names for pySCENIC
  "--cell_types",
  action = "extend",
  nargs = "+",
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
  args <- parser$parse_args(c(
    "--ame_csv", "/CBL_scRNAseq/results/human/Aldinger/20230224/enriched_motifs_nonhomol.csv",
    "--rss_csv", "/CBL_scRNAseq/results/pyscenic/20230402/aldinger_RL.rss.csv",
    "--cell_types", "19-NonHomol UBC", "20-NonHomol UBC",
    "--out_dir", "/CBL_scRNAseq/results/human/Aldinger/20230420/"
  ))
} else {
  args <- parser$parse_args()
}

message(sprintf("Saving files to %s", args$out_dir))
if (!dir.exists(args$out_dir)) {
  dir.create(args$out_dir, recursive = TRUE)
}

# ------------------------------------------------------------------------------
# import functions

source("/CBL_scRNAseq/software/utilities/pyscenic_regulons.R")
source("/CBL_scRNAseq/software/utilities/plot_venn_diagrams.R")

# ------------------------------------------------------------------------------
# create list of enriched TF motifs from the AME output

# import AME data
ame_motifs <- read_csv(args$ame_csv)

# parse TF names from `motif_id` column
ame_tfs <- ame_motifs$motif_id %>% 
  str_remove_all(pattern = "_.*")

message("Saving AME TFs to ame_tfs.txt")
write_lines(ame_tfs, file = file.path(args$out_dir, "ame_tfs.txt"))

# ------------------------------------------------------------------------------
# create list of regulons from pySCENIC

# import regulon specificity scores
rss <- read.csv(args$rss_csv, row.names = 1)

# get top 50% of regulons from non-homologous UBCs
# function from `pyscenic_regulons.R`
pyscenic_tfs <- map(
  .x = args$cell_types,
  \(.x) get_top_regulons(cell_type = .x, rss = rss)
)

if (length(args$cell_types > 1)) {
  # combine the TFs from all cell clusters
  pyscenic_tfs <- pyscenic_tfs %>% unlist() %>% unique()
}

# save pyscenic TFs to file
message("Saving pySCENIC TFs to ame_tfs.txt")
write_lines(pyscenic_tfs, file = file.path(args$out_dir, "pyscenic_tfs.txt"))

# ------------------------------------------------------------------------------
# make Venn diagram

venn_data <- list(ame_tfs, pyscenic_tfs) %>% 
  `names<-`(c("AME", "pySCENIC")) %>% 
  Venn() %>% 
  process_data()

# function from `plot_venn_diagrams.R`
plt <- venn_plot(venn_data) +
  # additional customizations
  scale_x_continuous(expand = expansion(mult = 0.5)) + 
  scale_y_continuous(expand = expansion(mult = c(0.25, 0.1))) + 
  scale_fill_distiller(palette = "Blues", direction = 1, guide = "none") + 
  theme_void()

ggsave(
  file = "ame_pyscenic_genes_venn.png",
  plot = plt,
  path = args$out_dir,
  width = 6,
  height = 4,
  units = "in",
  dpi = 600
)
