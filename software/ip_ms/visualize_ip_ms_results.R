# ==============================================================================
# Generate volcano plot of the BRCA1 and Kaiso interactors.
# ==============================================================================

# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(argparse)
library(tidyverse)
library(ggrepel)

# parse command line arguments
parser <- ArgumentParser()

parser$add_argument(
  # name of Excel file
  "--filename",
  default = NULL,
  required = TRUE
)
parser$add_argument(
  # name of Excel sheets to read in
  "--sheets",
  action = "extend",
  default = NULL,
  nargs = "+"
)
parser$add_argument(
  # output directory
  "--out_dir",
  default = NULL,
  required = TRUE
)

if (interactive()) {
  arg_list <- parser$parse_args(c(
    "--filename", "/isilon/CBL_scRNAseq-archived/data/ip_ms/BRCA1_Kaiso_interaction_results_20240102.xlsx",
    "--sheets", "BRCA1_SAINTID6835", "Kaiso_SAINTID6836",
    "--out_dir", "/CBL_scRNAseq/results/ip_ms/20240105/"
  ))
} else {
  arg_list <- parser$parse_args()
}

# load functions
source("/CBL_scRNAseq/software/ip_ms/load_ip_ms_data.R")
source("/CBL_scRNAseq/software/utilities/plotting.R")

# ------------------------------------------------------------------------------
# load data

dat <- map(
  .x = arg_list$sheets,
  .f = \(x) {
    load_ip_ms_data(
      filename = arg_list$filename,
      sheet = x
    )
  }
)
names(dat) <- arg_list$sheets

# ------------------------------------------------------------------------------
# volcano plot

# add labels for plotting
for (i in 1:length(dat)) {
  dat[[i]] <- mutate(
    dat[[i]],
    label = case_when(BFDR < 0.01 ~ PreyGene,
                      TRUE ~ NA_character_)
  )
}

# for BRCA1
.plt <- make_volcano(
  data = dat$BRCA1_SAINTID6835,
  log_fc = log2(FoldChange),
  log_pval = -log10(BFDR + 1),
  direction = NULL,
  gene_labels = label,
  log_pval_thresh = -log10(0.01 + 1)
) + 
  labs(title = "BRCA1", x = "log2(FC)", y = "-log10(p-val + 1)") + 
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.25)))
ggsave(
  filename = "brca1_volcano.png",
  plot = .plt,
  path = arg_list$out_dir,
  width = 12,
  height = 6,
  units = "in",
  dpi = 600
)

# for Kaiso
.plt <- make_volcano(
  data = dat$Kaiso_SAINTID6836,
  log_fc = log2(FoldChange),
  log_pval = -log10(BFDR + 1),
  direction = NULL,
  gene_labels = label,
  log_pval_thresh = -log10(0.01 + 1)
) + 
  labs(title = "Kaiso", x = "log2(FC)", y = "-log10(p-val + 1)") + 
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.2)))
ggsave(
  filename = "kaiso_volcano.png",
  plot = .plt,
  path = arg_list$out_dir,
  width = 12,
  height = 6,
  units = "in",
  dpi = 600
)
