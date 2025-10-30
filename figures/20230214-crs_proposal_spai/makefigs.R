# ==============================================================================
# These figures were generated for Shraddha's CRS grant proposal due 2023-02-14
# ==============================================================================

# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(tidyverse)
library(Seurat)
library(ggrepel)

# Seurat objects for import
srat_names <- c("aldinger", "integ_RL", "integ_full")
srat_rds <- c(
  # Aldinger human RL from Liam
  "/isilon/CBL_scRNAseq-archived/data/human/Aldinger/glutamatergic_dev_Liam.RDS",
  # CCA - RL only
  "/CBL_scRNAseq/results/integrated/vladoiu_liam_RL.rds",
  # CCA - full cerebellum
  "/CBL_scRNAseq/results/integrated/20230126/without_future/aldinger_vladoiu_cca.rds"
)

# import Seurat objects
all_srat <- purrr::map(
  .x = srat_rds,
  .f = readRDS
)
names(all_srat) <- srat_names
print(all_srat)

# import functions
source("/CBL_scRNAseq/software/mouse/Vladoiu/add_annotations.R")
source("/CBL_scRNAseq/software/utilities/cell_labelling.R")
source("/CBL_scRNAseq/software/utilities/cluster_barplot.R")
source("/CBL_scRNAseq/software/utilities/score_integration.R")

# ------------------------------------------------------------------------------
# plot number of human vs. mouse cells in each cluster

# make plot (`cluster_barplot.R`)
.plt <- cluster_barplot(
  all_srat$integ_RL,
  split.by = "species",
  group.by = "seurat_clusters",
  position = "stack",
  filter_data = "seurat_clusters %in% c(7, 19, 20)"
) + 
  labs(x = "Cell cluster", y = "Number of cells") + 
  scale_fill_manual(values = c("#027FA3", "#BFBFBF")) + 
  theme(legend.position = "none")

ggsave(
  filename = "num_cells_barplot.png",
  plot = .plt,
  width = 2.5,
  height = 3,
  units = "in",
  dpi = 600
)

# ------------------------------------------------------------------------------
# plot distance between cluster centroids

# add mouse cell types (`add_annotation.R`)
all_srat$integ_full <- label_vladoiu_cells(all_srat$integ_full) %>% 
  # pool cell types (`cell_labelling.R`)
  pool_cell_types(col_name = "common_cell_type")

# get embeddings and metadata
embeddings <- Embeddings(all_srat$integ_full, reduction = "pca")
metadata <- all_srat$integ_full@meta.data

# calculate distance between human and mouse centroids for each cell type
all_dist <- delta_centroids(
  embeddings = embeddings,
  metadata = metadata,
  by = "common_cell_type",
  filter_out = c("other/missing")
)

# fix order of integ_method and common_cell_type
all_dist <- all_dist %>% 
  # remove astrocytes, microglia, GCPs, GNs from graph
  filter(
    !(common_cell_type %in% c("astrocytes", "microglia", "GCPs", "GNs"))
  ) %>% 
  # reorder so that the OPC controls show up at the end
  mutate(
    common_cell_type = fct_relevel(common_cell_type, "RL", "GCPs", "GNs", "UBCs", "OPCs")
  )

# make plot
.plt <- ggplot(all_dist, aes(x = common_cell_type, y = delta)) + 
  geom_point(size = 3) + 
  labs(
    x = "Cell type",
    y = "Distance between human\nand mouse cluster centroids",
    colour = "Integration\nmethod",
    shape = "Integration\nmethod"
  ) + 
  scale_y_continuous(limits = c(0, NA)) + # `NA` uses the current min/max
  theme_light() + 
  theme(axis.text = element_text(colour = "black"))

ggsave(
  filename = "delta_centroids.png",
  plot = .plt,
  width = 3,
  height = 3,
  units = "in",
  dpi = 600
)

# ------------------------------------------------------------------------------
# volcano plot of differentially expressed genes in homologous vs non-homologous UBCs

de_genes <- read_csv("/CBL_scRNAseq/results/human/Aldinger/20230205/de_genes.csv")

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

# make plot
.plt <- ggplot(
  data = de_genes,
  aes(
    x = avg_log2FC,
    y = -log(p_val_adj, base = 10),
    label = gene_label
  )
) + 
  geom_point(aes(colour = diff_exp), size = 1) + 
  geom_text_repel() + 
  geom_vline(xintercept = c(logfc_threshold, -logfc_threshold), colour = "red", linewidth = 0.25) + 
  geom_hline(yintercept = -log(0.05, base = 10), colour = "red", linewidth = 0.25) + 
  labs(x = expression("log"[2]*"(fold change)"),
       y = expression("-log"[10]*"(adjusted p-value)")) + 
  scale_colour_manual(
    values = c("blue", "red"),
    breaks = c("down", "up"),
    labels = c("Non-homologous\nUBC up", "Homologous\nUBC up"),
    name = NULL
  ) + 
  guides(colour = guide_legend(override.aes = list(size = 2))) + 
  theme_light() + 
  theme(
    axis.text = element_text(colour = "black"),
    legend.position = "bottom",
    legend.box.spacing = unit(0.05, "in")
  )

ggsave(
  "volcano_plot.png",
  plot = .plt,
  width = 4,
  height = 4,
  units = "in",
  dpi = 600
)
