# ==============================================================================
# Get differentially expressed genes in our non-homologous UBC population
# ==============================================================================

# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(argparse)
library(tidyverse)
library(Seurat)
library(ggrepel)

# parse command line arguments
parser <- ArgumentParser()

parser$add_argument(
  # path to human-only Seurat object for differential expression
  "--srat_rds",
  default = NULL,
  required = TRUE
)
parser$add_argument(
  # path to integrated Seurat object to copy homologous/non-homologous UBC labels
  "--integ_rds",
  default = NULL,
  required = TRUE
)
parser$add_argument(
  # output directory
  "--out_dir",
  default = NULL,
  required = TRUE
)

args <- parser$parse_args()

if (is.null(args$srat_rds)) {
  stop("Argument for `srat_rds` is missing; please provide a Seurat file for differential expression analysis")
} else {
  srat_rds <- args$srat_rds
}

if (is.null(args$integ_rds)) {
  stop("Argument for `integ_rds` is missing; please provide an integrated Seurat object file to copy the labels from")
} else {
  integ_rds <- args$integ_rds
}

if (is.null(args$out_dir)) {
  stop("Argument for `out_dir` is missing; please provide an output directory")
} else {
  out_dir <- args$out_dir
  if (!dir.exists(out_dir)) {
    dir.create(out_dir)
  }
}

# ------------------------------------------------------------------------------
# import functions

source("/CBL_scRNAseq/software/utilities/cell_labelling.R")

# ------------------------------------------------------------------------------
# import data

message(sprintf("Reading in from: %s", srat_rds))
diffexp_srat <- readRDS(srat_rds)

message(sprintf("Reading in from: %s", integ_rds))
integ_srat <- readRDS(integ_rds)

message(sprintf("Saving output to: %s", out_dir))

# add cell type to integ_clusters
integ_srat <- label_integ_clusters(integ_srat, integ_method = "cca")

# copy the integrated cluster numbers to the human dataset, making sure that the
# cells have the correct cluster labelled
diffexp_srat$integ_clusters <- integ_srat@meta.data[row.names(diffexp_srat@meta.data), "integ_clusters"]

# ------------------------------------------------------------------------------
# compare gene expression in homologous vs non-homologous UBCs

diffexp_srat <- PrepSCTFindMarkers(diffexp_srat)
Idents(diffexp_srat) <- "integ_clusters"

de_genes <- lapply(
  X = list(c("7-Homol UBC", "19-NonHomol UBC", "20-NonHomol UBC")),
  FUN = function(X) {
    df <- FindMarkers(
      diffexp_srat,
      assay = "SCT",
      slot = "data",
      logfc.threshold = 0, # default 0.25
      test.use = "wilcox",
      min.pct = 0.1, # default 0.1
      ident.1 = X[1],
      ident.2 = X[-1],
      # group.by = "seurat_clusters",
      # subset.ident = "human", # take only cells in these clusters
      # recorrect_umi = FALSE # required if subsetting cells from original object
    ) %>% 
      mutate(gene = rownames(.),
             which_clusters = paste(X[1], 
                                    paste(X[-1], collapse = "+"), 
                                    sep = " vs. "),
             .before = 1) %>% 
      `rownames<-`(., NULL)
  }
) %>% 
  bind_rows(.)

# save to file
readr::write_csv(de_genes, file = file.path(out_dir, "de_genes.csv"))

# ------------------------------------------------------------------------------
# make volcano plot

# set threshold for differential expression
logfc_threshold <- 1

# add labels
de_genes <- mutate(
  de_genes,
  diff_exp = case_when(
    (avg_log2FC <= -logfc_threshold & p_val_adj < 0.05) ~ "down",
    (avg_log2FC >= logfc_threshold & p_val_adj < 0.05) ~ "up"
  ),
  gene_label = case_when(!is.na(diff_exp) ~ gene)
)

plt <- ggplot(
  data = de_genes[de_genes$which_clusters == "7-Homol UBC vs. 19-NonHomol UBC+20-NonHomol UBC", ],
  aes(
    x = avg_log2FC,
    y = -log(p_val_adj, base = 10),
    label = gene_label
  )
) + 
  geom_point(size = 1) + 
  geom_text_repel() + 
  geom_vline(xintercept = c(logfc_threshold, -logfc_threshold), colour = "red", size = 0.25) + 
  geom_hline(yintercept = -log(0.05, base = 10), colour = "red", size = 0.25) + 
  labs(x = expression("log"[2]*"(fold change)"),
       y = expression("-log"[10]*"(adjusted p-value)")) + 
  theme_minimal() + 
  theme(axis.text = element_text(colour = "black"))
ggsave("volcano_plot.png", plot = plt, path = out_dir, width = 8, height = 8, units = "in")

# ------------------------------------------------------------------------------

print(sessionInfo())
