# ==============================================================================
# These figures were generated for Ian's committee meeting #1 on 
# Monday, April 4, 2022
# ==============================================================================

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(networkD3)
library(tidyverse)
library(Seurat)

# import data
integ_srat <- readRDS("/CBL_scRNAseq/results/integrated/vladoiu_liam_integrated.rds")
vlad_srat <- readRDS("/CBL_scRNAseq/results/mouse/Vladoiu/merged_seurat.rds")
liam_srat <- readRDS("/isilon/CBL_scRNAseq/data/human/Aldinger/glutamatergic_dev_Liam.RDS")

# U of T colour palette
uoft_colours <- c("#1E3765", # navy blue
                  "#007FA3", # secondary blue
                  "#6D247A", # purple
                  "#DC4633", # red
                  "#00A189", # blue-green
                  "#6FC7EA", # sky blue
                  "#AB1368", # dark magenta
                  "#0D534D", # dark green
                  "#F1C500", # yellow
                  "#8DBF2E", # almost OICR green
                  "#000000", # black
                  "#D0D1C9") # light grey

# ==============================================================================
# make alluvial/sankey chart of orthologous genes
# ==============================================================================

# nodes
nodes <- data.frame(genes = c("mouse genes (18638)", 
                              "1:1 ortholog (14544)", 
                              "no orthologs found (3360)", 
                              "multiple orthologs (732)", 
                              "manually excluded (2)"))

# edges
links <- as.data.frame(matrix(c(0, 1, 14544, 
                                0, 2, 3360, 
                                0, 3, 732, 
                                0, 4, 2), 
                              byrow = TRUE, ncol = 3))
names(links) <- c("source", "target", "value")

# make Sankey network
orthologs <- sankeyNetwork(Links = links, 
                           Nodes = nodes, 
                           Source = "source", 
                           Target = "target", 
                           Value = "value", 
                           NodeID = "genes", 
                           fontSize = 72, 
                           fontFamily = "Helvetica", 
                           nodeWidth = 48, 
                           nodePadding = 96, 
                           height = 1080, 
                           width = 1920)

# save as html and then convert to png
orthologs_html <- "./orthologs.html"
saveNetwork(orthologs, file = orthologs_html)
# webshot::webshot(url = orthologs_html, file = "./orthologs.png", 
#                  vwidth = 1920, vheight = 1080)

# ==============================================================================
# Vladoiu and Aldinger QC, etc.
# ==============================================================================



# ==============================================================================
# UMAP of Vladoiu and Aldinger
# ==============================================================================

DimPlot(vlad_srat, reduction = "umap", label = TRUE, group.by = "orig.ident") + NoLegend()
ggsave(filename = "vlad_umap.png", width = 5, height = 5, units = "in", dpi = 600)

DimPlot(liam_srat, reduction = "umap", label = TRUE) + NoLegend()
ggsave(filename = "liam_umap.png", width = 5, height = 5, units = "in", dpi = 600)

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
  annotate("segment", x = 7.05, xend = 7.05, y = 4000, yend = 3600, size = 0.75, 
           arrow = arrow(length = unit(8, "bigpts"))) + 
  annotate("text", x = 7.1, y = 4400, label = "EOMES+\n(UBCs)") + 
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
# gene expression dot plot
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

# ==============================================================================
# gene expression feature plot of Purkinje markers, Aldinger data
# ==============================================================================

FeaturePlot(liam_srat, features = c("CALB1", "CA8", "RORA"), order = TRUE, ncol = 3)
ggsave("purkinje_markers_human.png", width = 12, height = 4, units = "in", dpi = 600)