# ==============================================================================
# These figures were generated for Ian's student seminar on 2023-03-07.
# ==============================================================================

# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(tidyverse)
library(Seurat)
library(slingshot)
library(monocle3)
library(ggrepel)

# Seurat objects for import
srat_names <- c("aldinger", "vladoiu", "cca_RL", "cca_full")
srat_rds <- c(
  # Aldinger human RL from Liam
  "/isilon/CBL_scRNAseq-archived/data/human/Aldinger/glutamatergic_dev_Liam.RDS",
  # Vladoiu mouse RL
  "/CBL_scRNAseq/results/mouse/Vladoiu/merged_seurat_RLonly.rds",
  # CCA RL only
  "/CBL_scRNAseq/results/integrated/vladoiu_liam_RL.rds",
  # CCA full human/mouse cerebellum
  "/CBL_scRNAseq/results/integrated/20230126/without_future/aldinger_vladoiu_cca.rds"
)

# import Seurat objects
all_srat <- purrr::map(
  .x = srat_rds,
  .f = readRDS
)
names(all_srat) <- srat_names
print(all_srat)

# import slingshot pseudotime result
# sling_pto <- readRDS(
#   "/CBL_scRNAseq/results/human/Aldinger/20221124/umap.NA.new_cell_type.RL-VZ/pto.rds"
# )

# import monocle pseudotime results
# monocle_cds <- readRDS("/CBL_scRNAseq/results/human/Aldinger/20221209/monocle_cds.rds")

# import functions
source("/CBL_scRNAseq/software/mouse/Vladoiu/add_annotations.R")
source("/CBL_scRNAseq/software/utilities/cell_labelling.R")
source("/CBL_scRNAseq/software/utilities/score_integration.R")
source("/CBL_scRNAseq/software/utilities/cluster_barplot.R")
# source("/CBL_scRNAseq/software/utilities/slingshot.R")

# colour palette for species
species_cols <- c(hcl(h = 15, c = 100, l = 65), "grey")

# ------------------------------------------------------------------------------
# UMAPs of full CCA integration

# label vladoiu cells(`add_annotations.R`)
all_srat$cca_full <- label_vladoiu_cells(all_srat$cca_full)

# pool human and mouse cell types (`cell_labelling.R`)
all_srat$cca_full <- pool_cell_types(all_srat$cca_full)

# cell names
human_cells <- rownames(all_srat$cca_full@meta.data)[all_srat$cca_full$species == "human"]
mouse_cells <- rownames(all_srat$cca_full@meta.data)[all_srat$cca_full$species == "mouse"]

.plt <- DimPlot(
  all_srat$cca_full,
  reduction = "umap",
  group.by = "seurat_clusters",
  label = TRUE,
  repel = TRUE,
  raster = FALSE
) + 
  NoLegend() + 
  labs(title = NULL)

ggsave(
  "umap_cca_clusters.png",
  plot = .plt,
  width = 5.5,
  height = 5,
  units = "in",
  dpi = 600
)

.plt <- DimPlot(
  all_srat$cca_full,
  group.by = "species",
  pt.size = 0.01,
  cells.highlight = human_cells,
  sizes.highlight = 0.01,
  label = FALSE,
  raster = FALSE
) + 
  scale_colour_manual(labels = c("mouse", "human"), values = rev(species_cols)) +
  labs(title = NULL)

ggsave(
  "umap_cca_species.png",
  plot = .plt,
  width = 6.25,
  height = 5,
  units = "in",
  dpi = 600
)

umap_human_cells <- DimPlot(
  all_srat$cca_full,
  cells = human_cells,
  reduction = "umap",
  group.by = "figure_clusters",
  label = TRUE,
  label.size = 2.5,
  repel = TRUE,
  raster = FALSE
) + 
  NoLegend() + 
  ggtitle("human cell types")

umap_mouse_cells <- DimPlot(
  all_srat$cca_full,
  cells = mouse_cells,
  cols = c(
    # RL/UBC lineage
    scales::brewer_pal(palette = "Reds")(4),
    # granule cell (progenitors)
    scales::brewer_pal(palette = "Blues")(3),
    # other glutamatergic
    scales::brewer_pal(palette = "YlOrBr")(3),
    # stem cells
    scales::brewer_pal(palette = "BuPu")(5),
    # GABAergic lineage
    scales::brewer_pal(palette = "Greens")(7),
    # glial/non-neuronal cells
    scales::brewer_pal(palette = "PuRd")(9),
    # NA
    "#DFDFDF"
  ),
  reduction = "umap",
  group.by = "mouse_cell_type",
  label = TRUE,
  label.size = 2.5,
  repel = TRUE,
  raster = FALSE
) + 
  NoLegend() + 
  labs(title = "mouse cell types")

ggsave(
  filename = "umap_celltypes.png",
  plot = umap_human_cells + umap_mouse_cells,
  width = 11,
  height = 5,
  units = "in",
  dpi = 600
)

# OPCs
human_opcs <- rownames(all_srat$cca_full@meta.data)[all_srat$cca_full$figure_clusters %in% c("11-OPC", "12-Committed OPC")]
mouse_opcs <- rownames(all_srat$cca_full@meta.data)[all_srat$cca_full$mouse_cell_type %in% c("Oligodendrocyte precursor cells")]

umap_human_opcs <- DimPlot(
  all_srat$cca_full,
  cells = human_cells,
  reduction = "umap",
  raster = FALSE,
  cells.highlight = human_opcs,
  # cols.highlight = "#1E3765",
  sizes.highlight = 0.1
) + 
  NoLegend() + CenterTitle() + 
  labs(title = "human OPCs") + 
  scale_colour_manual(values = c(hcl(15, 50, 80), "#1E3765"))

umap_mouse_opcs <- DimPlot(
  all_srat$cca_full,
  cells = mouse_cells,
  reduction = "umap",
  raster = FALSE,
  cells.highlight = mouse_opcs,
  cols.highlight = "#1E3765",
  sizes.highlight = 0.1
) + 
  NoLegend() + CenterTitle() + 
  labs(title = "mouse OPCs")

ggsave(
  filename = "umap_oligodendrocytes.png",
  plot = umap_human_opcs + umap_mouse_opcs,
  width = 11,
  height = 5,
  units = "in",
  dpi = 600
)

# ------------------------------------------------------------------------------
# plot distance between cluster centroids

# get embeddings and metadata
embeddings <- Embeddings(all_srat$cca_full, reduction = "pca")
metadata <- all_srat$cca_full@meta.data

# filter out cell types (don't calculate distance for these)
common_cell_types <- unique(all_srat$cca_full$common_cell_type)
filter_out <- common_cell_types[!(common_cell_types %in% c("RL", "GCPs", "GNs", "UBCs", "OPCs"))]

# calculate distance between human and mouse centroids for each cell type
all_dist <- delta_centroids(
  embeddings = embeddings,
  metadata = metadata,
  by = "common_cell_type",
  filter_out = filter_out
)

# fix order of integ_method and common_cell_type
all_dist <- all_dist %>% 
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
  scale_y_continuous(limits = c(10, NA)) + # `NA` uses the current min/max
  theme_light() + 
  theme(axis.text = element_text(colour = "black"),
        axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(
  filename = "delta_centroids.png",
  plot = .plt,
  width = 2.5,
  height = 4.5,
  units = "in",
  dpi = 1200
)

# ------------------------------------------------------------------------------

all_srat$cca_full@meta.data <- mutate(
  all_srat$cca_full@meta.data,
  common_cell_type = case_when(
    common_cell_type %in% c("RL", "GCPs", "GNs", "UBCs", "OPCs") ~ common_cell_type,
    TRUE ~ "other"
  ) %>% fct_relevel("RL", "GCPs", "GNs", "UBCs", "OPCs", "other")
)

.plt <- DimPlot(
  all_srat$cca_full,
  group.by = "common_cell_type",
  cols = c(scales::hue_pal()(4), "#1E3765", "grey"),
  label = TRUE,
  repel = TRUE,
  raster = FALSE
) + 
  NoLegend() + 
  labs(title = NULL)

ggsave(
  filename = "umap_common_cells.png",
  plot = .plt,
  width = 5.5,
  height = 5,
  units = "in",
  dpi = 600
)

# ------------------------------------------------------------------------------

# import differential gene expression results
de_genes <- read_csv("/CBL_scRNAseq/results/human/Aldinger/20230205/de_genes.csv")

# set threshold for differential expression
logfc_threshold <- 1

# add labels
# human-specific UBC genes get swapped for the plot ()
de_genes <- mutate(
  de_genes,
  avg_log2FC = -1 * avg_log2FC,
  diff_exp = case_when(
    (avg_log2FC < 0 & p_val_adj < 0.05) ~ "down",
    (avg_log2FC > 0 & p_val_adj < 0.05) ~ "up"
  ),
  gene_label = case_when(abs(avg_log2FC) > 1 ~ gene)
)

# make volcano plot
.plt <- ggplot(
  data = de_genes[de_genes$which_clusters == "7-Homol UBC vs. 19-NonHomol UBC+20-NonHomol UBC", ],
  aes(
    x = avg_log2FC,
    y = -log(p_val_adj, base = 10),
    label = gene_label,
    colour = diff_exp
  )
) + 
  geom_point(size = 1) + 
  geom_text_repel() + 
  # geom_vline(xintercept = c(logfc_threshold, -logfc_threshold), colour = "red", size = 0.25) + 
  geom_hline(yintercept = -log(0.05, base = 10), colour = "black", linewidth = 0.5, linetype = "dashed") + 
  labs(x = expression("log"[2]*"(fold change)"),
       y = expression("-log"[10]*"(adjusted p-value)")) + 
  scale_colour_manual(values = c("blue", "red")) + 
  theme_classic() + 
  theme(
    axis.text = element_text(colour = "black"),
    legend.position = "none"
  )

ggsave(
  "volcano_plot.png",
  plot = .plt,
  width = 4,
  height = 4,
  units = "in",
  dpi = 600
)

# ------------------------------------------------------------------------------
# barplot of mouse vs human cells in each cluster

# relevel human and mouse so human shows on bottom of plot
all_srat$cca_RL$species <- fct_relevel(all_srat$cca_RL$species, "mouse", "human")

.plt <- cluster_barplot(
  all_srat$cca_RL,
  split.by = "species",
  width = 0.6,
  filter_data = "seurat_clusters %in% c(7, 19, 20)"
) + 
  xlab("cluster") + 
  scale_fill_manual(labels = c("mouse", "human"), values = rev(species_cols)) + 
  theme_classic() + 
  theme(
    legend.position = c(0.85, 0.9),
    axis.text = element_text(colour = "black"),
    axis.ticks = element_line(colour = "black")
  )

ggsave(
  "species_per_cluster.png",
  plot = .plt,
  width = 3.5,
  height = 4,
  units = "in",
  dpi = 600
)

# ------------------------------------------------------------------------------

print(sessionInfo())
