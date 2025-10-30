# ==============================================================================
# Get shared pySCENIC regulons and MB-mutated genes.
# ==============================================================================

# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(argparse)
library(tidyverse)
library(ggVennDiagram)
library(patchwork)

# parse command line arguments
parser <- ArgumentParser()

parser$add_argument(
  # path to CSV of RSS scores
  "--rss_csv",
  default = NULL,
  required = TRUE
)
parser$add_argument(
  # number of regulons to compare
  "--num_regulons",
  default = NULL,
  required = TRUE,
  type = "integer"
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
    "--rss_csv", "/CBL_scRNAseq/results/pyscenic/20230320/aldinger_RL.rss.csv",
    "--num_regulons", 8,
    "--out_dir", "/CBL_scRNAseq/results/pyscenic/testing/"
  ))
} else {
  args <- parser$parse_args()
}

message(sprintf("Saving files to %s", args$out_dir))
if (!dir.exists(args$out_dir)) {
  dir.create(args$out_dir, recursive = TRUE)
}

# ------------------------------------------------------------------------------
# load top pySCENIC regulons using the regulon specificity scores

# import data
rss <- read.csv(args$rss_csv, row.names = 1)

cell_types <- rownames(rss)

# top regulons for each cluster/cell type
top_regulons <- lapply(
  X = as.list(cell_types),
  FUN = function(x, rss, num_regulons) {
    # get the RSS scores
    scores <- data.matrix(rss)[x, ]
    
    # sort scores then get the TF associated with the score
    top_regs <- scores %>% 
      sort(decreasing = TRUE) %>% 
      names() %>% 
      head(n = num_regulons)
    
    return(top_regs)
  },
  rss = rss, # RSS scores
  num_regulons = args$num_regulons # top n regulons
)

names(top_regulons) <- cell_types

# ------------------------------------------------------------------------------
# load MB-mutated genes (see lab notebook 20230223)

source("/CBL_scRNAseq/software/human/Hendrikse/load_data.R")

# function from `load_data.R`
mut_genes <- load_mut_genes(
  mut_xlsx = "/isilon/CBL_scRNAseq-archived/data/human/Hendrikse/supp_tab_5.xlsx", # default
  rename_cols = TRUE # default
)

# filter for significantly mutated genes
sig_mut_genes <- mut_genes %>% 
  mutate(
    g4_mutated = case_when(
      qval_oncodrive_g4 < 0.05 | qval_mutsig_g4 < 0.05 ~ TRUE,
      TRUE ~ FALSE
    ),
    g3_mutated = case_when(
      qval_oncodrive_g3 < 0.05 | qval_mutsig_g3 < 0.05 ~ TRUE,
      TRUE ~ FALSE
    ),
    mutated_in = case_when(
      g4_mutated & g3_mutated ~ "both",
      g4_mutated ~ "G4",
      g3_mutated ~ "G3",
      TRUE ~ NA_character_
    )
  ) %>% 
  select(gene, Freq_whole_cohort, freq_g4, freq_g3, mutated_in) %>% 
  filter(!is.na(mutated_in))

# additional genes from Figure 1 of Hendrikse et al. (2022)
additional_genes <- c(
  "PRDM6",
  "CBFA2T3",
  "KDM6A",
  "RUNX1T1",
  "CBFA2T2",
  "KDM2B",
  "KBTBD4",
  "SMARCA4",
  "KMT2D",
  "KMT2C",
  "SETD2",
  "SRCAP",
  "ZMYM3",
  "EHMT1",
  "CREBBP",
  "EHMT2",
  "KAT6A",
  "KDM6B",
  "KMT2A",
  "CHD7",
  "CHD9",
  "CHD8",
  "CHD3",
  "CHD4",
  "CHD6",
  "CHD2",
  "CHD5",
  "GFI1B",
  "GFI1",
  "OTX2",
  "ZIC1",
  "TBR1",
  "MYC",
  "MYCN",
  "PVT1",
  "CTDNEP1",
  "CDK6",
  "PKD1",
  "STARD9",
  "ASPM",
  "HTT",
  "ANKRD17",
  "EIF2AK4",
  "TSC2",
  "ATM",
  "FANCA",
  "FANCL",
  "PALB2",
  "FANCI",
  "BRCA2",
  "FANCB",
  "FANCG",
  "FANCM",
  "ZFHX2",
  "ZFHX3",
  "ZFHX4",
  "IKBKAP",
  "ELP2",
  "ELP4",
  "DST"
)

all_mut_genes <- union(
  sig_mut_genes$gene,
  additional_genes
)

# ------------------------------------------------------------------------------
# Venn diagram for common TFs in cerebellar development and mutated genes in MB.

venn_plts <- map(
  .x = cell_types,
  .f = function(x, top_regulons, mb_genes) {
    # subset top TFs for this specific cell type
    tfs <- top_regulons[[x]]
    
    # combine TFs and mutated genes into list and prep data for plotting
    venn_data <- list(tfs, mb_genes) %>% 
      `names<-`(c(x, "mb_mut_genes")) %>% 
      Venn() %>% 
      process_data()
    
    # get intersection of TFs and mutated genes
    common_genes <- dplyr::intersect(tfs, mb_genes)
    
    # make the plot
    .plt <- ggplot() + 
      
      # the set areas
      geom_sf(aes(fill = count), data = venn_region(venn_data)) + 
      
      # the set edges
      geom_sf(color = "black", size = 2, data = venn_setedge(venn_data), show.legend = FALSE) + 
      
      # name of the sets
      geom_sf_text(
        aes(label = name),
        data = venn_setlabel(venn_data),
        size = 6,
        fontface = "bold"
      ) +
      
      # counts
      geom_sf_text(aes(label = count), data = venn_region(venn_data), size = 5) + 
      
      # label intersecting genes
      geom_segment(
        aes(x = 500, y = 400, xend = 500, yend = 200),
        color = "red",
        linewidth = 1
      ) + 
      geom_text(
        aes(x = 500, y = 180, label = paste(common_genes, collapse = "\n")),
        data = venn_region(venn_data),
        size = 5,
        vjust = 1
      ) + 
      
      # additional customizations
      scale_x_continuous(expand = expansion(mult = 0.5)) + 
      scale_y_continuous(expand = expansion(mult = c(0.25, 0.1))) + 
      scale_fill_distiller(palette = "Blues", direction = 1, guide = "none") + 
      theme_void()
    
    return(.plt)
  },
  top_regulons = top_regulons,
  mb_genes = all_mut_genes
)

plt_dim <- length(cell_types) %>% sqrt() %>% ceiling()
venn_plts <- wrap_plots(
  venn_plts,
  ncol = plt_dim,
  widths = 2,
  heights = 1
)

ggsave(
  filename = "tfs_mb_mut_shared_genes_venn.png",
  plot = venn_plts,
  path = args$out_dir,
  width = plt_dim * 6,
  height = plt_dim * 3,
  units = "in",
  dpi = 300
)

# ------------------------------------------------------------------------------

message("\nSESSION INFO\n")
print(sessionInfo())

