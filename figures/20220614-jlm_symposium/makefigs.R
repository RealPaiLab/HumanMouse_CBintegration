# ------------------------------------------------------------------------------
# These figures were generated for Ian's virtual poster presentation at the
# MBP/JLM Symposium on Tuesday, June 14, 2022.
# ------------------------------------------------------------------------------

# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(patchwork)
library(Seurat)
library(tidyverse)

# import the integrated Vladoiu + Aldinger RL-only Seurat object
integ_srat <- readRDS("/CBL_scRNAseq/results/integrated/vladoiu_liam_RL.rds")

# ------------------------------------------------------------------------------
# UMAPs
# ------------------------------------------------------------------------------

# colour by clusters
umap_cluster <- DimPlot(integ_srat, label = TRUE, repel = TRUE) + NoLegend()

# colour by species
species_cols <- c("#007FA3", "grey")
umap_species <- DimPlot(integ_srat, label = FALSE, group.by = "species", pt.size = 0.01,
                        cells.highlight = rownames(integ_srat@meta.data)[integ_srat$species == "human"],
                        sizes.highlight = 0.01) +
  scale_colour_manual(labels = c("mouse", "human"),
                      values = rev(species_cols)) + 
  labs(title = NULL) +
  theme(legend.direction = "horizontal",
        legend.position = c(0.5, 1),
        legend.justification = c("center", "bottom"))

# ------------------------------------------------------------------------------
# EOMES expression on UMAP
# ------------------------------------------------------------------------------

eomes_expr <- FeaturePlot(integ_srat, features = "EOMES", 
                          min.cutoff = "q10", max.cutoff = "q90", 
                          order = TRUE)

# ------------------------------------------------------------------------------
# bar plot showing proportion of human/mouse cells in each cluster
# ------------------------------------------------------------------------------

num_cells <- table(integ_srat$species, integ_srat$seurat_clusters) %>% 
  as.data.frame(.)
colnames(num_cells) <- c("species", "cluster", "freq")

# plot absolute number of human/mouse cells in the UBC clusters (8, 17, 18)
num_cells_species <- ggplot(num_cells[num_cells$cluster %in% c(8, 17, 18), ], 
                            aes(x = cluster, y = freq, fill = fct_rev(species))) + 
  geom_col(width = 0.5) + 
  labs(y = "number") + 
  theme_light() + 
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) + 
  scale_fill_manual(values = rev(species_cols), name = NULL) + 
  theme_classic() + 
  theme(legend.direction = "horizontal",
        legend.position = c(0.5, 1.05),
        legend.justification = c("center", "bottom"),
        axis.text = element_text(colour = "black"), 
        axis.ticks = element_line(colour = "black"))

# ------------------------------------------------------------------------------
# combine all four plots
# ------------------------------------------------------------------------------

# set background fill
bg_fill <- function() {return(element_rect(fill = alpha("white", alpha = 0), colour = NA))}

# combine plots
fig1 <- (umap_cluster | umap_species) / (eomes_expr | num_cells_species) + 
  plot_annotation(tag_levels = "A") & 
  theme(panel.background = bg_fill(), 
        plot.background = bg_fill(), 
        legend.background = bg_fill(), 
        plot.tag = element_text(face = "bold"), 
        plot.tag.position = "topleft", 
        text = element_text(size = 15))

ggsave("fig1.png", plot = fig1, width = 20, height = 17, units = "cm")
ggsave("fig1.svg", plot = fig1, width = 20, height = 17, units = "cm", fix_text_size = FALSE)
