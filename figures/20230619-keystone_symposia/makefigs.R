# ==============================================================================
# These figures were generated for the Keystone Symposia conference in June
# 2023. It contains figures for both the poster and the short talk.
# ==============================================================================

# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(tidyverse)
library(patchwork)
library(ggrepel)
library(Seurat)
library(slingshot)
library(monocle3)

# Seurat objects for import
srat_names <- c("aldinger", "vladoiu", "cca_RL", "cca_full")
srat_rds <- c(
  # Aldinger human RL from Liam
  "/isilon/CBL_scRNAseq-archived/data/human/Aldinger/glutamatergic_dev_Liam.RDS",
  # Vladoiu mouse full cerebellum
  "/CBL_scRNAseq/results/mouse/Vladoiu/merged_seurat.rds",
  # CCA RL only
  "/CBL_scRNAseq/results/integrated/vladoiu_liam_RL.rds",
  # CCA full cerebellum
  "/CBL_scRNAseq/results/integrated/20230126/without_future/aldinger_vladoiu_cca.rds"
)

# import Seurat objects
all_srat <- purrr::map(
  .x = srat_rds,
  .f = readRDS
)
names(all_srat) <- srat_names
print(all_srat)

# import monocle pseudotime results
monocle_cds <- readRDS("/CBL_scRNAseq/results/human/Aldinger/20221209/monocle_cds.rds")

# import functions
source("/CBL_scRNAseq/software/mouse/Vladoiu/add_annotations.R")
source("/CBL_scRNAseq/software/utilities/cell_labelling.R")
# source("/CBL_scRNAseq/software/utilities/score_integration.R")
source("/CBL_scRNAseq/software/utilities/cluster_barplot.R")
# source("/CBL_scRNAseq/software/utilities/slingshot.R")

# colour palette for species
# species_cols <- c(hcl(h = 15, c = 100, l = 65), "#95D5D1")
species_cols <- c("#65BC45", "#95D5D1")

# convert to factor
all_srat$cca_RL$species <- fct_relevel(all_srat$cca_RL$species, "mouse", "human")
# all_srat$cca_full$species <- fct_relevel(all_srat$cca_full$species, "mouse", "human")

# adding labels to fully integrated data (functions from add_annotations.R and
# cell_labelling.R)
# generates a metadata column called `common_cell_type`
all_srat$cca_full <- label_vladoiu_cells(all_srat$cca_full) %>% 
  pool_cell_types(.)

# group cell types into OPC, RL-derived, or other
all_srat$cca_full@meta.data <- mutate(
  all_srat$cca_full@meta.data,
  common_cell_type = case_when(
    common_cell_type == "OPCs" ~ common_cell_type,
    common_cell_type %in% c("RL", "GCPs", "GNs", "UBCs") ~ "RL-derived",
    TRUE ~ "other"
  ) %>% fct_relevel("OPCs", "RL-derived", "other")
)

################################################################################
# for poster

out_dir <- file.path(".", "poster")

# ------------------------------------------------------------------------------
# human/mouse UMAP pre-integration

# mouse RL-derived cells
mouse_rl <- purrr::map(
  .x = c(0, 1, 3, 6, 9, 11, 13, 16, 21, 22, 34),
  .f = \(x) {
    rownames(all_srat$vladoiu@meta.data)[all_srat$vladoiu$seurat_clusters %in% x]
  }
)

# UMAP of human RL-derived cells
.plt1 <- DimPlot(
  all_srat$aldinger,
  reduction = "umap",
  group.by = "new_cell_type",
  cells.highlight = rownames(all_srat$aldinger@meta.data),
  sizes.highlight = 0.01,
  label = FALSE,
) + 
  NoLegend() +
  scale_colour_manual(labels = "human", values = species_cols[1]) + 
  labs(title = NULL) + 
  theme(axis.text = element_blank(), axis.ticks = element_blank())

# UMAP of mouse RL-derived cells
.plt2 <- DimPlot(
  all_srat$vladoiu,
  reduction = "umap",
  group.by = "seurat_clusters",
  cells.highlight = mouse_rl,
  cols.highlight = species_cols[2],
  sizes.highlight = 0.01,
  label = FALSE,
) + 
  NoLegend() + 
  labs(title = NULL) + 
  theme(axis.text = element_blank(), axis.ticks = element_blank())

# UMAP of integrated cells
.plt3 <- DimPlot(
  all_srat$cca_RL,
  reduction = "umap",
  group.by = "species",
  cols = species_cols,
  label = FALSE,
  repel = TRUE
) + 
  NoLegend() + 
  labs(title = "CCA integration") + 
  theme(axis.text = element_blank(), axis.ticks = element_blank())

layout <- c(
  area(t = 1, l = 1, r = 3),
  area(t = 1, l = 5, r = 7),
  area(t = 2, l = 3, r = 5)
)
.plt <- .plt1 + .plt2 + .plt3 + 
  plot_layout(design = layout) + 
  plot_annotation(tag_levels = 'A')

ggsave(
  filename = "1-integration.png",
  plot = .plt,
  path = out_dir,
  width = 8,
  height = 6,
  units = "in",
  dpi = 1200
)

# ------------------------------------------------------------------------------
# integration figures

# cell names
human_cells <- rownames(all_srat$cca_RL@meta.data)[all_srat$cca_RL$species == "human"]
mouse_cells <- rownames(all_srat$cca_RL@meta.data)[all_srat$cca_RL$species == "mouse"]
ubcs <- purrr::map(
  .x = c(7, 19, 20),
  .f = \(x) {
    rownames(all_srat$cca_RL@meta.data)[all_srat$cca_RL$seurat_clusters %in% x]
  }
)

# CCA integration showing clustering
.plt1 <- DimPlot(
  all_srat$cca_RL,
  reduction = "umap",
  group.by = "seurat_clusters",
  cells.highlight = ubcs,
  cols.highlight = rev(RColorBrewer::brewer.pal(n = 3, name = "Dark2")),
  sizes.highlight = 0.01,
  label = TRUE,
  repel = TRUE
) + 
  NoLegend() + 
  labs(title = "Three UBC clusters") + 
  theme(axis.text = element_blank(), axis.ticks = element_blank())

# EOMES on UMAP
.plt2 <- FeaturePlot(
  all_srat$cca_RL,
  features = "EOMES",
  min.cutoff = "q10",
  max.cutoff = "q90",
  order = TRUE
) + 
  labs(title = "*EOMES* (UBC marker)") + 
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    plot.title = ggtext::element_markdown()
  )

# bar plot - species per cluster
.plt3 <- cluster_barplot(
  all_srat$cca_RL,
  split.by = "species",
  width = 0.6,
  filter_data = "seurat_clusters %in% c(7, 19, 20)"
) + 
  labs(x = "clusters", fill = NULL) + 
  scale_fill_manual(labels = c("mouse", "human"), values = rev(species_cols)) + 
  theme_classic() + 
  theme(
    legend.position = c(0.7, 0.9),
    axis.text = element_text(colour = "black", size = 11),
    axis.ticks = element_line(colour = "black")
  )

# UMAP of full cerebellum integration with OPCs for control
human_cells <- rownames(all_srat$cca_full@meta.data)[all_srat$cca_full$species == "human"]
.plt4 <- DimPlot(
  all_srat$cca_full,
  group.by = "species",
  pt.size = 0.01,
  cells.highlight = human_cells,
  sizes.highlight = 0.01,
  # cols = c("#C33C54", "#DEA3A0", "grey"),
  label = FALSE,
  raster = FALSE
) + 
  NoLegend() + 
  scale_colour_manual(labels = c("mouse", "human"), values = rev(species_cols)) +
  labs(title = "Whole cerebellum\nwith OPC control") + 
  theme(axis.text = element_blank(), axis.ticks = element_blank())

# combine plots
layout <- "AABBCDD"
.plt <- .plt2 + .plt1 + .plt3 + .plt4 + 
  plot_layout(design = layout) + 
  plot_annotation(tag_levels = "A") & 
  theme(
    plot.tag = element_text(face = "bold", size = 14, vjust = 1),
    plot.tag.position = c(0, 1),
    plot.title = element_blank(),
    axis.title = element_text(size = 14)
  )

ggsave(
  filename = "2-human_spec_ubcs.png",
  plot = .plt,
  path = out_dir,
  width = 12,
  height = 3,
  units = "in",
  dpi = 1200
)

# ------------------------------------------------------------------------------
# differential expression, enriched TFs

# import differential gene expression results
de_genes <- read_csv("/CBL_scRNAseq/results/human/Aldinger/20230205/de_genes.csv")

# set threshold for differential expression
logfc_threshold <- 1

# add labels
# human-specific UBC genes get swapped for the plot
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
  geom_text_repel(size = 3) + 
  # geom_vline(xintercept = c(logfc_threshold, -logfc_threshold), colour = "red", size = 0.25) + 
  geom_hline(yintercept = -log(0.05, base = 10), colour = "black", linewidth = 0.5, linetype = "dashed") + 
  labs(x = expression("log"[2]*"(fold change)"),
       y = expression("-log"[10]*"(adjusted p-value)")) + 
  scale_color_brewer(palette = "Set1", direction = -1) + 
  theme_classic() + 
  theme(
    axis.text = element_text(colour = "black"),
    legend.position = "none",
    plot.tag = element_text(face = "bold", vjust = 1)
  )

ggsave(
  "3-ubc_diff_exp.png",
  plot = .plt,
  path = out_dir,
  width = 2.5,
  height = 2.5,
  units = "in",
  dpi = 1200
)



# import list of AME TFs
enr_motifs <- read_csv("/CBL_scRNAseq/results/human/Aldinger/20230224/enriched_motifs_nonhomol.csv")

# subset GFI1, GFI1B, PRDM6, MYCN, ZIC1
enr_motifs_sub <- enr_motifs %>% 
  mutate(
    motif_id = str_remove_all(string = motif_id, pattern = "_.*")
  ) %>% 
  filter(
    motif_id %in% c("GFI1", "GFI1B", "PRDM6", "MYCN", "ZIC1")
  ) %>% 
  select(motif_id, adj.pvalue)

.plt <- memes::plot_ame_heatmap(enr_motifs_sub) + 
  # add numbers
  geom_text(
    aes(x = motif_id, y = "All Regions", label = round(-log10(adj.pvalue), 1)),
    data = enr_motifs_sub
  ) + 
  labs(
    x = "Enrichment of known MB drivers",
    tag = "C",
    fill = "-log<sub>10</sub>(adj.<br>p-value)"
  ) + 
  scale_x_discrete(expand = expansion(mult = 0), limits = rev) + 
  scale_y_discrete(expand = expansion(mult = 0)) + 
  coord_flip() + 
  theme(
    axis.ticks.x = element_blank(),
    axis.text.x = element_blank(),
    legend.title = ggtext::element_markdown(),
    plot.tag = element_text(face = "bold")
  )

ggsave(
  "4-key_enr_motifs.png",
  plot = .plt,
  path = out_dir,
  width = 2.75,
  height = 3.5,
  units = "in",
  dpi = 1200
)

# ------------------------------------------------------------------------------
# monocle pseudotime

.plt <- plot_cells(
  monocle_cds,
  color_cells_by = "integ_clusters",
  label_groups_by_cluster = FALSE,
  group_label_size = 4
) + 
  theme_void() + 
  theme(legend.position = "none")

ggsave(
  "5-pseudotime.png",
  plot = .plt,
  path = out_dir,
  width = 4,
  height = 3,
  units = "in",
  dpi = 1200
)

################################################################################
# for talk

out_dir <- file.path(".", "talk")

# ------------------------------------------------------------------------------
# full cerebellar integration

# coloured by species
human_cells <- rownames(all_srat$cca_full@meta.data)[all_srat$cca_full$species == "human"]
.plt <- DimPlot(
  all_srat$cca_full,
  group.by = "species",
  pt.size = 0.01,
  cells.highlight = human_cells,
  sizes.highlight = 0.01,
  label = FALSE,
  raster = FALSE
) + 
  NoLegend() + 
  scale_colour_manual(labels = c("mouse", "human"), values = rev(species_cols)) +
  labs(title = NULL) + 
  theme(axis.text = element_blank(), axis.ticks = element_blank())

ggsave(
  filename = "full_integ_umap_species.png",
  plot = .plt,
  path = out_dir,
  width = 3.25,
  height = 3,
  units = "in",
  dpi = 600
)

# OPCs highlighted
.plt <- DimPlot(
  all_srat$cca_full,
  group.by = "common_cell_type",
  cols = c("#1E3765", scales::hue_pal()(1), "grey"),
  label = TRUE,
  repel = FALSE,
  raster = FALSE
) + 
  NoLegend() + 
  labs(title = NULL) + 
  theme(axis.text = element_blank(), axis.ticks = element_blank())

ggsave(
  filename = "full_integ_umap_opcs.png",
  plot = .plt,
  path = out_dir,
  width = 3.25,
  height = 3,
  units = "in",
  dpi = 600
)



# ------------------------------------------------------------------------------
# integration figures

# mouse RL-derived cells
mouse_rl <- purrr::map(
  .x = c(0, 1, 3, 6, 9, 11, 13, 16, 21, 22, 34),
  .f = \(x) {
    rownames(all_srat$vladoiu@meta.data)[all_srat$vladoiu$seurat_clusters %in% x]
  }
)

# UMAP of human RL-derived cells
.plt <- DimPlot(
  all_srat$aldinger,
  reduction = "umap",
  group.by = "new_cell_type",
  cells.highlight = rownames(all_srat$aldinger@meta.data),
  sizes.highlight = 0.01,
  label = FALSE,
) + 
  NoLegend() +
  scale_colour_manual(labels = "human", values = species_cols[1]) + 
  labs(title = NULL) + 
  theme(axis.text = element_blank(), axis.ticks = element_blank())

ggsave(
  filename = "human_umap.png",
  plot = .plt,
  path = out_dir,
  width = 2.5,
  height = 2.5,
  units = "in",
  dpi = 600
)

# UMAP of mouse RL-derived cells
.plt <- DimPlot(
  all_srat$vladoiu,
  reduction = "umap",
  group.by = "seurat_clusters",
  cells.highlight = mouse_rl,
  cols.highlight = species_cols[2],
  sizes.highlight = 0.01,
  label = FALSE,
) + 
  NoLegend() + 
  labs(title = NULL) + 
  theme(axis.text = element_blank(), axis.ticks = element_blank())

ggsave(
  filename = "mouse_umap.png",
  plot = .plt,
  path = out_dir,
  width = 2.5,
  height = 2.5,
  units = "in",
  dpi = 600
)

# UMAP of integrated cells
human_cells <- rownames(all_srat$cca_RL@meta.data)[all_srat$cca_RL$species == "human"]
.plt <- DimPlot(
  all_srat$cca_RL,
  group.by = "species",
  pt.size = 0.01,
  cells.highlight = human_cells,
  sizes.highlight = 0.01,
  label = FALSE,
  raster = FALSE
) + 
  NoLegend() + 
  scale_colour_manual(labels = c("mouse", "human"), values = rev(species_cols)) +
  labs(title = NULL) + 
  theme(axis.text = element_blank(), axis.ticks = element_blank())

ggsave(
  filename = "integ_umap_species.png",
  plot = .plt,
  path = out_dir,
  width = 2.75,
  height = 2.5,
  units = "in",
  dpi = 600
)

# ------------------------------------------------------------------------------
# EOMES expression

.plt <- FeaturePlot(
  all_srat$cca_RL,
  features = "EOMES",
  order = TRUE,
  min.cutoff = "q10",
  max.cutoff = "q90"
) + 
  labs(title = "*EOMES* (UBC marker)") + 
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    plot.title = ggtext::element_markdown()
  )

ggsave(
  filename = "integ_umap_eomes.png",
  plot = .plt,
  path = out_dir,
  width = 3.75,
  height = 3,
  units = "in",
  dpi = 600
)

# ------------------------------------------------------------------------------
# integration Seurat clusters

.plt <- DimPlot(
  all_srat$cca_RL,
  reduction = "umap",
  group.by = "seurat_clusters",
  label = TRUE,
  repel = TRUE
) + 
  NoLegend() + 
  labs(title = NULL) + 
  theme(axis.text = element_blank(), axis.ticks = element_blank())

ggsave(
  filename = "integ_umap_clusters.png",
  plot = .plt,
  path = out_dir,
  width = 3.25,
  height = 3,
  units = "in",
  dpi = 600
)

# ------------------------------------------------------------------------------
# bar plot - species per cluster

.plt <- cluster_barplot(
  all_srat$cca_RL,
  split.by = "species",
  width = 0.6,
  filter_data = "seurat_clusters %in% c(7, 19, 20)"
) + 
  labs(x = "clusters", fill = NULL) + 
  scale_fill_manual(labels = c("mouse", "human"), values = rev(species_cols)) + 
  theme_classic() + 
  theme(
    legend.position = c(0.65, 0.9),
    axis.text = element_text(colour = "black"),
    axis.ticks = element_line(colour = "black")
  )

ggsave(
  "human_enriched_clusters_bar.png",
  plot = .plt,
  path = out_dir,
  width = 2.25,
  height = 3,
  units = "in",
  dpi = 600
)

# ------------------------------------------------------------------------------
# differential expression

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
  geom_text_repel(size = 3) + 
  # geom_vline(xintercept = c(logfc_threshold, -logfc_threshold), colour = "red", size = 0.25) + 
  geom_hline(yintercept = -log(0.05, base = 10), colour = "black", linewidth = 0.5, linetype = "dashed") + 
  labs(x = expression("log"[2]*"(fold change)"),
       y = expression("-log"[10]*"(adjusted p-value)")) + 
  scale_color_brewer(palette = "Set1", direction = -1) + 
  theme_classic() + 
  theme(
    axis.text = element_text(colour = "black"),
    axis.ticks = element_line(colour = "black"),
    legend.position = "none",
  )

ggsave(
  "volcano_human_ubcs.png",
  plot = .plt,
  path = out_dir,
  width = 3,
  height = 3,
  units = "in",
  dpi = 600
)

# ------------------------------------------------------------------------------
# pseudotime

.plt <- plot_cells(
  monocle_cds,
  color_cells_by = "pseudotime"
) + 
  theme_void()

ggsave(
  "pseudotime.png",
  plot = .plt,
  path = out_dir,
  width = 4,
  height = 3,
  units = "in",
  dpi = 1200
)

# ------------------------------------------------------------------------------

print(sessionInfo())
