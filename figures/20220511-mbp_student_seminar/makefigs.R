# ==============================================================================
# These figures were generated for Ian's MBP Student Seminar on Wednesday,
# May 11, 2022
# ==============================================================================

# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(tidyverse)
library(Seurat)

# import data; this is the integration of Vladoiu (mouse RL-derived cells only)
# and Liam/Aldinger (human RL-derived cells only) as opposed to previous figures
# which contained the entire Vladoiu dataset
integ_srat <- readRDS("/CBL_scRNAseq/results/integrated/vladoiu_liam_RL.rds")

# import Vladoiu data
mouse_srat <- readRDS("/CBL_scRNAseq/results/mouse/Vladoiu/merged_seurat.rds")

# ==============================================================================
# Vladoiu clusters selected for integration (RL-derived cells)
# ==============================================================================

# get RL-derived cells (see notes from 2022-04-21 for cluster info)
rl_cells <- WhichCells(mouse_srat, idents = c(1, 3, 5:8, 11:12, 16, 19, 33))

# UMAP highlighting the RL-derived cells
rl_cells_umap <- DimPlot(mouse_srat, 
                         cells.highlight = rl_cells, 
                         sizes.highlight = 0.01) + 
  scale_colour_manual(labels = c("non-RL-derived", "RL-derived"), values = c("grey", "#AB1368"))
ggsave("mouse_umap_RLcells.png", plot = rl_cells_umap,
       width = 6.5, height = 5, units = "in", dpi = 600)

# UMAP of all mouse cells
mouse_umap <- DimPlot(mouse_srat, label = TRUE, repel = TRUE) + NoLegend()
ggsave("mouse_umap_clusters.png", plot = mouse_umap,
       width = 5.5, height = 5, units = "in", dpi = 600)

# both UMAPs on same plot
mouse_umap + rl_cells_umap
ggsave("mouse_umaps.png", width = 12, height = 5, units = "in", dpi = 600)

# ==============================================================================
# UMAP of integrated data
# ==============================================================================

# set two-tone palette
species_cols <- c(hcl(h = 15, c = 100, l = 65), "grey")

# UMAP by species
species_umap <- DimPlot(integ_srat, label = FALSE, group.by = "species", pt.size = 0.01, 
                        cells.highlight = rownames(integ_srat@meta.data)[integ_srat$species == "human"], 
                        sizes.highlight = 0.01) + 
  scale_color_manual(labels = c("mouse", "human"), values = rev(species_cols)) + 
  labs(title = NULL)
ggsave("integrated_umap_species.png", plot = species_umap, 
       width = 6.5, height = 5, units = "in", dpi = 600)

# UMAP by clusters
cluster_umap <- DimPlot(integ_srat, label = TRUE, repel = TRUE) + NoLegend()
ggsave("integrated_umap_clusters.png", plot = cluster_umap, 
       width = 5.5, height = 5, units = "in", dpi = 600)

# both UMAPs on same plot
cluster_umap + species_umap
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
  # annotate("segment", x = c(9.05, 18.05, 19.05), xend = c(9.05, 18.05, 19.05),
  #          y = 4000, yend = 3600, size = 0.75, 
  #          arrow = arrow(length = unit(8, "bigpts"))) + 
  # annotate("text", x = 7.1, y = 4400, label = "EOMES+\n(UBCs)") + 
  labs(y = "number") + 
  theme_light() + 
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) + 
  scale_fill_manual(values = rev(species_cols), name = "species") + 
  theme_classic() + 
  theme(legend.position = c(0.9, 0.5), 
        axis.text = element_text(colour = "black"), 
        axis.text.x = element_text(angle = 45, hjust = 1), 
        axis.ticks = element_line(colour = "black"))
ggsave("species_bar.png", width = 9, height = 5, units = "in", dpi = 600)

# ==============================================================================
# gene expression dot plot and UMAP
# ==============================================================================

# major mouse gene markers from Carter et al. 2018 and Vladoiu et al. 2019
mouse_markers <- c(
  "Atoh1", # upper rhombic lip, glutamatergic
  "Wls", # "interior" rhombic lip (Yeung et al. 2014)
  "Ptf1a", # ventricular zone, GABAergic
  "Hes5", "Id1", "Msx3", "Nes", "Sox2", # progenitor/neural stem cells
  "Msx1", # roof plate
  "Lmx1a", # roof plate and unipolar brush cells
  "Eomes", "Calb2",  # unipolar brush cells
  "Pax6", # granule neuron progenitors
  "Meis2", "Lhx2", # glutamatergic cerebellar nuclei/nuclear transitory neurons
  "Tbr1", # glutamatergic/excitatory cerebellar nuclei
  "Calb1", "Car8", "Rora", # Purkinje cells
  "Pax2", "Lbx1", # GABAergic interneurons
  "Sox10", "Olig1" # oligodendrocytes
)

# major human gene markers
m2h_genes <- read.csv("/CBL_scRNAseq/results/integrated/vladoiu_orth_genes.csv")
m2h_genes <-  m2h_genes[duplicated(m2h_genes$MGI.symbol) == FALSE
                        & duplicated(m2h_genes$HGNC.symbol) == FALSE
                        & m2h_genes$MGI.symbol != "Pisd", ]
human_markers <- m2h_genes[match(mouse_markers, m2h_genes$MGI.symbol), ]$HGNC.symbol %>% 
  .[!is.na(.)]

# make dot plot
DotPlot(integ_srat, features = human_markers, assay = "integrated") + 
  scale_x_discrete(limits = rev) + 
  coord_flip()
ggsave("gene_expression_dotplot.png", width = 12, height = 6, units = "in", dpi = 600)

# UBC marker expression on UMAP
ubc_genes <- c("EOMES", "LMX1A")
for (gene in ubc_genes) {
  FeaturePlot(integ_srat, features = gene, min.cutoff = "q10", max.cutoff = "q90", order = TRUE)
  ggsave(paste0(gene, ".png"), width = 5, height = 4, units = "in", dpi = 600)
}

# ==============================================================================
# Liam's annotated cells in each cluster
# ==============================================================================

human_annot <- table(integ_srat$new_cell_type, integ_srat$seurat_clusters) %>% 
  as.data.frame(.)
colnames(human_annot) <- c("annot", "cluster", "freq")

# plot number of annotated human cells in each cluster
ggplot(human_annot, aes(x = cluster, y = freq, fill = annot)) + 
  geom_col() + 
  labs(y = "number") + 
  theme_light() + 
  scale_y_continuous(expand = expansion(mult = c(0, 0.05)))
ggsave("human_cells_per_cluster.png", width = 12, height = 5, units = "in", dpi = 600)
