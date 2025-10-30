# ==============================================================================
# These figures were generated for Ian's committee meeting #3 on Monday,
# September 25, 2023.
# ==============================================================================

# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(tidyverse)
library(patchwork)
library(ggalluvial)
library(Seurat)
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

# convert to factor
all_srat$cca_RL$species <- fct_relevel(all_srat$cca_RL$species, "mouse", "human")
all_srat$cca_full$species <- fct_relevel(all_srat$cca_full$species, "mouse", "human")

# import monocle pseudotime results
rl_cds <- readRDS("/CBL_scRNAseq/results/human/Aldinger/20230718/monocle_cds.rds")
ubc_cds <- readRDS("/CBL_scRNAseq/results/human/Aldinger/20230719/monocle_cds_ubc_lineage.rds")

# colour palette for species
species_cols <- c(hcl(h = 15, c = 100, l = 65), "grey")

# import functions
source("/CBL_scRNAseq/software/mouse/Vladoiu/add_annotations.R")
source("/CBL_scRNAseq/software/utilities/cell_labelling.R")
source("/CBL_scRNAseq/software/utilities/plotting.R")

# ------------------------------------------------------------------------------
# alluvial diagram of orthologs

orth_levels <- c("1:1 ortholog",
                 "no orthologs found",
                 "multiple orthologs",
                 "manually excluded")
orthologs <- tribble(
  ~mouse, ~genes, ~human, ~freq,
  "mouse", "1:1 ortholog", "human", 14544,
  "mouse", "no orthologs found", NA, 3360,
  "mouse", "multiple orthologs", NA, 732,
  "mouse", "manually excluded", NA, 2,
  NA, "no orthologs found", "human", 1e3,
  NA, "multiple orthologs", "human", 1e2,
  NA, "manually excluded", "human", 1e1,
) %>% 
  to_lodes_form(axes = c("mouse", "genes", "human")) %>% 
  mutate(
    orth_status = case_when(
      alluvium == 1 ~ "1:1 ortholog",
      alluvium %in% c(2, 5) ~ "no orthologs found",
      alluvium %in% c(3, 6) ~ "multiple orthologs",
      alluvium %in% c(4, 7) ~ "manually excluded",
    ),
    orth_status = factor(orth_status, levels = orth_levels)
  )

.plt <- ggplot(
  data = orthologs,
  mapping = aes(x = x, y = freq, alluvium = alluvium, stratum = stratum)
) + 
  geom_alluvium(aes(fill = orth_status), na.rm = TRUE) + 
  geom_stratum(na.rm = TRUE)

# ------------------------------------------------------------------------------
# bar chart of human/mouse cells from each cluster

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
    legend.position = "right",
    axis.text = element_text(colour = "black", size = 11),
    axis.ticks = element_line(colour = "black")
  )

ggsave(
  "human_enriched_clusters_bar.png",
  plot = .plt,
  width = 5,
  height = 3,
  units = "in",
  dpi = 600
)

# ------------------------------------------------------------------------------
# UMAP of full cerebellar integration with OPCs as controls

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
  width = 3.25,
  height = 3,
  units = "in",
  dpi = 600
)

# OPCs highlighted
.plt <- DimPlot(
  all_srat$cca_full,
  group.by = "common_cell_type",
  cols = c("#1E3765", "#FDCDAC", "#B3E2CD"),
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
  width = 3.25,
  height = 3,
  units = "in",
  dpi = 600
)

# ------------------------------------------------------------------------------
# pseudotime plots

# >>> for full RL >>>

# by pseudotime
.plt1 <- plot_cells(
  rl_cds,
  color_cells_by = "pseudotime",
  trajectory_graph_segment_size = 2,
  cell_size = 0.75,
  label_cell_groups = FALSE
) + 
  theme_void()

# by integrated clusters
.plt2 <- plot_cells(
  rl_cds,
  color_cells_by = "integ_clusters",
  trajectory_graph_segment_size = 2,
  cell_size = 0.75,
  label_cell_groups = FALSE
) + 
  scale_colour_manual(
    values = DiscretePalette(
      n = length(levels(colData(rl_cds)$integ_clusters)),
      palette = "polychrome"
    )
  ) + 
  theme_void() + 
  guides(color = guide_legend(
    ncol = 1,
    override.aes = list(size = 2),
    title = NULL
  ))

ggsave(
  "pseudotime_full_rl.png",
  plot = (.plt1 + .plt2),
  width = 15,
  height = 5,
  units = "in"
)

# <<<

# >>> for UBC lineage only >>>

# by pseudotime
.plt1 <- plot_cells(
  ubc_cds,
  color_cells_by = "pseudotime",
  trajectory_graph_segment_size = 2,
  cell_size = 0.75,
  label_cell_groups = FALSE
) + 
  theme_void()

# by integrated clusters
.plt2 <- plot_cells(
  ubc_cds,
  color_cells_by = "integ_clusters",
  trajectory_graph_segment_size = 2,
  cell_size = 0.75,
  label_cell_groups = FALSE
) + 
  scale_colour_manual(
    values = DiscretePalette(
      n = length(levels(colData(ubc_cds)$integ_clusters)),
      palette = "polychrome"
    )
  ) + 
  theme_void() + 
  guides(color = guide_legend(
    ncol = 1,
    override.aes = list(size = 2),
    title = NULL
  ))

ggsave(
  "pseudotime_ubc_lineage.png",
  plot = (.plt1 + .plt2),
  width = 15,
  height = 5,
  units = "in"
)

# <<<


# ------------------------------------------------------------------------------
# species-specific differential expression

# >>> rl-svz >>>

rlsvz <- read_csv("/CBL_scRNAseq/results/species_dge/20230829/02-rl_svz.csv")

.plt <- make_volcano(
  data = rlsvz,
  log_fc = logFC,
  log_pval = -log10(FDR),
  log_pval_thresh = -log10(0.05)
) + 
  labs(y = "-log10 FDR") + 
  theme_classic() + 
  theme(
    axis.text = element_text(colour = "black"),
    axis.ticks = element_line(colour = "black")
  )

ggsave(
  "rl_svz_volcano.png",
  plot = .plt,
  width = 5,
  height = 5,
  units = "in"
)

# <<<

# >>> rl-vz >>>

rlvz <- read_csv("/CBL_scRNAseq/results/species_dge/20230829/01-rl_vz.csv")

.plt <- make_volcano(
  data = rlvz,
  log_fc = logFC,
  log_pval = -log10(FDR),
  log_pval_thresh = -log10(0.05)
) + 
  labs(y = "-log10 FDR") + 
  theme_classic() + 
  theme(
    axis.text = element_text(colour = "black"),
    axis.ticks = element_line(colour = "black")
  )

ggsave(
  "rl_vz_volcano.png",
  plot = .plt,
  width = 5,
  height = 5,
  units = "in"
)

# <<<
