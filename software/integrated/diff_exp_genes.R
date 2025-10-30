# ==============================================================================
# Get differentially expressed genes
# ==============================================================================

# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(dplyr)
library(ggplot2)
library(ggrepel)
library(Seurat)

# set data and output directories
root_dir <- "CBL_scRNAseq"
out_dir <- file.path(paste0("/", root_dir), "results/integrated")
date_dir <- file.path(out_dir, "20220715")

# ------------------------------------------------------------------------------
# load data

ubc_srat <- readRDS(file.path("", root_dir, "results/human/Aldinger/UBC_seurat.rds")) %>% 
  PrepSCTFindMarkers(.)

# ------------------------------------------------------------------------------
# compare human cells in clusters:
# 8 vs. 17 and 18
# 8 vs. 17
# 8 vs. 18
# 17 vs. 18

Idents(ubc_srat) <- "integ_clusters"
de_genes <- lapply(
  X = list(c(8, 17, 18),
           c(8, 17),
           c(8, 18),
           c(17, 18)),
  FUN = function(X) {
    df <- FindMarkers(
      ubc_srat,
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
readr::write_csv(de_genes, file = file.path(date_dir, "de_genes.csv"))

# ------------------------------------------------------------------------------
# make volcano plot of 8 vs. 17 and 18

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

ggplot(data = de_genes[de_genes$which_clusters == "8 vs. 17+18", ],
       aes(x = avg_log2FC, y = -log(p_val_adj, base = 10),
           label = gene_label)) + 
  geom_point(size = 1) + 
  geom_text_repel() + 
  geom_vline(xintercept = c(logfc_threshold, -logfc_threshold), colour = "red", size = 0.25) + 
  geom_hline(yintercept = -log(0.05, base = 10), colour = "red", size = 0.25) + 
  labs(x = expression("log"[2]*"(fold change)"),
       y = expression("-log"[10]*"(adjusted p-value)")) + 
  theme_minimal() + 
  theme(axis.text = element_text(colour = "black"))
ggsave("volcano_plot.pdf", path = date_dir, width = 6, height = 5, units = "in")

# ------------------------------------------------------------------------------
# save feature plots of the top differentially expressed genes

source(file.path("", root_dir, "software/utilities/check_feature_expression.R"))

# get only features from 8 vs. 17+18
features <- de_genes$gene_label[!is.na(de_genes$gene_label) & 
                                         de_genes$which_clusters == "8 vs. 17+18"]

# generate and save feature plots showing only human cells
dir.create(file.path(date_dir, "feature_plots"))
save_feature_plots(
  object = ubc_srat,
  features = features,
  plot_args = list(
    order = TRUE,
    min.cutoff = "q10",
    max.cutoff = "q90",
    label = TRUE,
    repel = TRUE
  ),
  save_args = list(
    path = file.path(date_dir, "feature_plots"),
    width = 5,
    height = 4,
    units = "in"
  )
)

# generate and save violin plots
dir.create(file.path(date_dir, "vln_plots"))
save_vln_plots(
  object = ubc_srat,
  features = features,
  plot_args = list(
    assay = "SCT"
  ),
  save_args = list(
    path = file.path(date_dir, "vln_plots"),
    width = 8,
    height = 5,
    units = "in"
  )
)
