# ==============================================================================
# Visualize and explore integrated data
# ==============================================================================

# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
# setwd("/CBL_scRNAseq/software/")

library(tidyverse)
library(Seurat)

# set data and output directories
root_dir <- "CBL_scRNAseq"
human_data_dir <- file.path("/isilon", root_dir, "data/human")
mouse_data_dir <- file.path("/isilon", root_dir, "data/mouse")

integ_out_dir <- file.path(paste0("/", root_dir), "results/integrated")

integ_date_dir <- file.path(integ_out_dir, format(Sys.Date(), "%Y%m%d"))

if (!dir.exists(integ_date_dir)) {
  dir.create(integ_date_dir)
}

# ==============================================================================
# import data
# ==============================================================================

integ_srat <- readRDS(file.path(integ_out_dir, "vladoiu_liam_integrated.rds"))

# ==============================================================================
# how many human/mouse cells are in each cluster?
# ==============================================================================

num_cells <- table(integ_srat$species, integ_srat$seurat_clusters) %>% 
  as.data.frame(.)
colnames(num_cells) <- c("species", "cluster", "freq")

# plot absolute number of human/mouse cells in each cluster
ggplot(num_cells, aes(x = cluster, y = freq, fill = species)) + 
  geom_col() + 
  labs(y = "number") + 
  theme_light() + 
  scale_y_continuous(expand = expansion(mult = c(0, 0.05)))
ggsave("num_cells_per_cluster.pdf", path = integ_date_dir, 
       width = 10, height = 5, units = "in")

# plot proportion of human/mouse cells in each cluster
ggplot(num_cells, aes(x = cluster, y = freq, fill = species)) + 
  geom_col(position = "fill") + 
  labs(y = "proportion") + 
  theme_light() + 
  scale_y_continuous(expand = expansion(mult = c(0, 0)))
ggsave("prop_cells_per_cluster.pdf", path = integ_date_dir,
       width = 10, height = 5, units = "in")

# ==============================================================================
# what marker genes are expressed in each cluster?
# ==============================================================================

# major mouse gene markers from Carter et al. 2018 and Vladoiu et al 2019
genes <- c(
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

# convert mouse genes to human orthologs
m2h_genes <- read_csv(file.path(integ_out_dir, "hgnc_mgi_orth_genes.csv"))
genes <- sapply(
  X = genes, 
  FUN = function(X) {
    if (X %in% m2h_genes$MGI.symbol) {
      # if an orthologous gene exists, then set the name to the new gene
      X <- m2h_genes$HGNC.symbol[m2h_genes$MGI.symbol == X]
    } else {
      # if no orthologous gene was found, keep the old (original) gene name
      X
    }
  }
)

# manually change mouse Msx3 to human VENTX then remove names
genes[names(genes) == "Msx3"] <- "VENTX"
genes <- unname(genes)

# add other human markers (Aldinger et al., supplementary table S8)
genes <- c(
  genes,
  "OTX2", # human RL_VZ marker
  "MKI67", # general proliferation marker
)

# make dot plot
DotPlot(integ_srat, features = genes, assay = "integrated") + scale_x_discrete(limits = rev) + coord_flip()
ggsave("gene_expression_dotplot.pdf", path = integ_date_dir, 
       width = 11.5, height = 8.5, units = "in")

# ==============================================================================
# based on Liam's annotations, what human cells are in each cluster?
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
ggsave("human_cells_per_cluster.pdf", path = integ_date_dir, 
       width = 10, height = 5, units = "in")

# ==============================================================================
# based on annotations from Vladoiu paper, what mouse cells are in each cluster?
# ==============================================================================

# change mouse cell types to factor and rename
source(file.path("", root_dir, "software/mouse/Vladoiu/add_annotations.R"))
integ_srat <- label_vladoiu_cells(integ_srat)

mouse_annot <- table(integ_srat$cell_type, integ_srat$seurat_clusters) %>% 
  as.data.frame(.)
colnames(mouse_annot) <- c("annot", "cluster", "freq")

# plot number of annotated mouse cells in each cluster
ggplot(mouse_annot, aes(x = cluster, y = freq, fill = annot)) + 
  geom_col() + 
  labs(y = "number") + 
  theme_light() + 
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) + 
  # use dittoSeq colour palette: 
  # https://github.com/dtm2451/dittoSeq/blob/master/R/dittoColors.R
  scale_fill_manual(values = c("#E69F00", "#56B4E9", "#009E73", "#F0E442",
                               "#0072B2", "#D55E00", "#CC79A7", "#666666",
                               "#AD7700", "#1C91D4", "#007756", "#D5C711",
                               "#005685", "#A04700", "#B14380", "#4D4D4D",
                               "#FFBE2D", "#80C7EF", "#00F6B3", "#F4EB71",
                               "#06A5FF", "#FF8320", "#D99BBD", "#8C8C8C",
                               "#FFCB57", "#9AD2F2", "#2CFFC6", "#F6EF8E",
                               "#38B7FF", "#FF9B4D", "#E0AFCA", "#A3A3A3",
                               "#8A5F00", "#1674A9", "#005F45", "#AA9F0D",
                               "#00446B", "#803800", "#8D3666", "#3D3D3D"))
ggsave("mouse_cells_per_cluster.pdf", path = integ_date_dir, 
       width = 20, height = 8, units = "in")
