# ==============================================================================
# These figures were generated for Shraddha's talk at the OICR Informatics
# Retreat on Friday, March 25, 2022.
# ==============================================================================

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(tidyverse)
library(Seurat)

# import data
integ_srat <- readRDS("/CBL_scRNAseq/results/integrated/vladoiu_liam_integrated.rds")

# set two-tone palette
species_cols <- c(hcl(h = 15, c = 100, l = 65), "grey")

# ==============================================================================
# make UMAP
# ==============================================================================

# UMAP by species
species_umap <- DimPlot(integ_srat, label = FALSE, group.by = "species", pt.size = 0.01, 
                        cells.highlight = rownames(integ_srat@meta.data)[integ_srat$species == "human"], 
                        sizes.highlight = 0.01) + 
  scale_color_manual(labels = c("mouse", "human"), values = rev(species_cols)) + 
  labs(title = NULL)

# UMAP by clusters
cluster_umap <- DimPlot(integ_srat, label = TRUE, repel = TRUE) + NoLegend()

species_umap + cluster_umap
ggsave("integrated_umaps.png", width = 12, height = 5, units = "in", dpi = 600)

# ==============================================================================
# make bar chart of human/mouse cells from each cluster
# ==============================================================================

num_cells <- table(integ_srat$species, integ_srat$seurat_clusters) %>% 
  as.data.frame(.)
colnames(num_cells) <- c("species", "cluster", "freq")

# plot absolute number of human/mouse cells in each cluster
ggplot(num_cells, aes(x = cluster, y = freq, fill = fct_rev(species))) + 
  geom_col() + 
  annotate("segment", x = 7.05, xend = 7.05, y = 4000, yend = 3600, size = 0.75, 
           arrow = arrow(length = unit(8, "bigpts"))) + 
  annotate("text", x = 7.1, y = 4400, label = "EOMES+\n(UBCs)") + 
  annotate("text", x = c(2, 24, 27), y = c(4300, 1400, 800), label = "*", size = 6) + 
  labs(y = "number") + 
  theme_light() + 
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) + 
  scale_fill_manual(values = rev(species_cols), name = "species") + 
  theme_classic() + 
  theme(legend.position = c(0.9, 0.5), 
        axis.text = element_text(colour = "black"), 
        axis.text.x = element_text(angle = 45, hjust = 1), 
        axis.ticks = element_line(colour = "black"))
ggsave("species_bar.png", width = 8, height = 4, units = "in", dpi = 600)
