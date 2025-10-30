# ==============================================================================
# These figures were generated for Shraddha's GRC poster in May 2023.
# ==============================================================================

# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(tidyverse)
library(Seurat)
library(monocle3)
library(ggrepel)
library(patchwork)

# Seurat objects for import
srat_names <- c("aldinger", "vladoiu", "cca_RL")
srat_rds <- c(
  # Aldinger human RL from Liam
  "/isilon/CBL_scRNAseq-archived/data/human/Aldinger/glutamatergic_dev_Liam.RDS",
  # Vladoiu mouse full cerebellum
  "/CBL_scRNAseq/results/mouse/Vladoiu/merged_seurat.rds",
  # CCA RL only
  "/CBL_scRNAseq/results/integrated/vladoiu_liam_RL.rds"
)

# import Seurat objects
all_srat <- purrr::map(
  .x = srat_rds,
  .f = readRDS
)
names(all_srat) <- srat_names
print(all_srat)

# import functions
# source("/CBL_scRNAseq/software/mouse/Vladoiu/add_annotations.R")
# source("/CBL_scRNAseq/software/utilities/cell_labelling.R")
# source("/CBL_scRNAseq/software/utilities/score_integration.R")
source("/CBL_scRNAseq/software/utilities/cluster_barplot.R")
# source("/CBL_scRNAseq/software/utilities/slingshot.R")

# colour palette for species
species_cols <- c(hcl(h = 15, c = 100, l = 65), "#A472BB")

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
  label = TRUE,
  repel = TRUE
) + 
  NoLegend() + 
  labs(title = NULL, tag = "A") + 
  theme(axis.text = element_blank(), axis.ticks = element_blank())

# UMAP of mouse RL-derived cells
.plt2 <- DimPlot(
  all_srat$vladoiu,
  reduction = "umap",
  group.by = "seurat_clusters",
  cells.highlight = mouse_rl,
  cols.highlight = scales::hue_pal()(length(mouse_rl)),
  sizes.highlight = 0.01,
  label = TRUE,
  repel = TRUE
) + 
  NoLegend() + 
  labs(title = NULL, tag = "B") + 
  theme(axis.text = element_blank(), axis.ticks = element_blank())

.plt <- .plt1 + .plt2 + 
  plot_layout(heights = c(1, 1))

ggsave(
  filename = "1-integration_a.png",
  plot = .plt,
  width = 3,
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
  labs(title = "Seurat clusters", tag = "D") + 
  theme(axis.text = element_blank(), axis.ticks = element_blank())

# EOMES on UMAP
.plt2 <- FeaturePlot(
  all_srat$cca_RL,
  features = "EOMES",
  min.cutoff = "q10",
  max.cutoff = "q90",
  order = TRUE
) + 
  labs(title = expression(bolditalic("EOMES")), tag = "C") + 
  theme(axis.text = element_blank(), axis.ticks = element_blank())

# bar plot - species per cluster
all_srat$cca_RL$species <- fct_relevel(all_srat$cca_RL$species, "mouse", "human")
.plt3 <- cluster_barplot(
  all_srat$cca_RL,
  split.by = "species",
  width = 0.6,
  filter_data = "seurat_clusters %in% c(7, 19, 20)"
) + 
  labs(x = "clusters", tag = "E", fill = NULL) + 
  scale_fill_manual(labels = c("mouse", "human"), values = rev(species_cols)) + 
  theme_classic() + 
  theme(
    legend.position = c(0.75, 0.85),
    axis.text = element_text(colour = "black"),
    axis.ticks = element_line(colour = "black")
  )

# combine plots
layout <- c(
  area(t = 1, b = 2, l = 1, r = 2),
  area(t = 3, b = 4, l = 1, r = 2),
  area(t = 1, b = 4, l = 3, r = 3)
)
.plt <- .plt2 + .plt1 + .plt3 + 
  plot_layout(design = layout) + 
  theme(plot.tag = element_text(face = "bold", size = 16, vjust = 0.5),
        plot.tag.position = c(0, 1))

ggsave(
  filename = "1-integration_b.png",
  plot = .plt,
  width = 6,
  height = 6,
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
  geom_text_repel() + 
  # geom_vline(xintercept = c(logfc_threshold, -logfc_threshold), colour = "red", size = 0.25) + 
  geom_hline(yintercept = -log(0.05, base = 10), colour = "black", linewidth = 0.5, linetype = "dashed") + 
  labs(x = expression("log"[2]*"(fold change)"),
       y = expression("-log"[10]*"(adjusted p-value)"),
       tag = "A") + 
  scale_colour_manual(values = rev(species_cols)) + 
  theme_classic() + 
  theme(
    axis.text = element_text(colour = "black"),
    legend.position = "none",
    plot.tag = element_text(face = "bold")
  )

ggsave(
  "2-ubc_diff_exp.png",
  plot = .plt,
  width = 3,
  height = 3,
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
  labs(x = "Enrichment of known MB drivers", tag = "C") + 
  scale_x_discrete(expand = expansion(mult = 0), limits = rev) + 
  scale_y_discrete(expand = expansion(mult = 0)) +
  coord_flip() + 
  theme(
    axis.ticks.x = element_blank(),
    axis.text.x = element_blank(),
    legend.title = element_blank(),
    plot.tag = element_text(face = "bold")
  )

ggsave(
  "3-key_enr_motifs.png",
  plot = .plt,
  width = 2.5,
  height = 3.5,
  units = "in",
  dpi = 1200
)

# ------------------------------------------------------------------------------
# monocle pseudotime

# import monocle pseudotime results
monocle_cds <- readRDS("/CBL_scRNAseq/results/human/Aldinger/20221209/monocle_cds.rds")

.plt <- plot_cells(
  monocle_cds,
  color_cells_by = "pseudotime",
  label_cell_groups = FALSE
) + 
  theme_void()
ggsave(
  "4-pseudotime.png",
  plot = .plt,
  width = 4,
  height = 3,
  units = "in",
  dpi = 1200
)

# ------------------------------------------------------------------------------

print(sessionInfo())
