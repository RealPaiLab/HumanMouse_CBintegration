# This script takes a Seurat object and generates basic visualizations,
# including gene expression and dot plots. This should be sourced from another
# file with the follwing variables already defined:
#   srat: a Seurat object
#   out_dir: the output directory
#   date_dir: the output directory with the date

library(Seurat)
library(tidyverse)


# ==============================================================================
# dot plot of genes 
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
m2h_genes <- read_csv("/CBL_scRNAseq/results/integrated/hgnc_mgi_orth_genes.csv")
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
  "MKI67" # general proliferation marker
)

# make dot plot
gene_expression_dotplot <-
  DotPlot(srat, features = genes, assay = "integrated") +
  scale_x_discrete(limits = rev) +
  coord_flip()
ggsave(filename = "gene_expression_dotplot.pdf",
       path = date_dir,
       width = 10, height = 7.5, units = "in")

cat(sprintf("saved dotplots to %s\n", date_dir))


# ==============================================================================
# bar chart of species in each cluster
# ==============================================================================

num_cells <- table(srat$species, srat$seurat_clusters) %>% 
  as.data.frame(.)
colnames(num_cells) <- c("species", "cluster", "freq")

# plot absolute number of human/mouse cells in each cluster
ggplot(num_cells, aes(x = cluster, y = freq, fill = species)) + 
  geom_col() + 
  labs(y = "number") + 
  theme_light() + 
  scale_y_continuous(expand = expansion(mult = c(0, 0.05)))
ggsave("num_cells_per_cluster.pdf", path = date_dir, 
       width = 10, height = 5, units = "in")

# plot proportion of human/mouse cells in each cluster
ggplot(num_cells, aes(x = cluster, y = freq, fill = species)) + 
  geom_col(position = "fill") + 
  labs(y = "proportion") + 
  theme_light() + 
  scale_y_continuous(expand = expansion(mult = c(0, 0)))
ggsave("prop_cells_per_cluster.pdf", path = date_dir,
       width = 10, height = 5, units = "in")

cat(sprintf("saved bar charts to %s\n", date_dir))
