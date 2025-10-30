# ==============================================================================
# Visualize and explore integrated data
# ==============================================================================

# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(argparse)
library(tidyverse)
library(Seurat)

# parse command line arguments
parser <- ArgumentParser()

parser$add_argument(
  # path to Seurat RDS object to import
  "--srat_rds",
  default = NULL,
  required = TRUE
)
parser$add_argument(
  # output directory
  "--out_dir",
  default = NULL,
  required = TRUE
)
parser$add_argument(
  # CSV file containing orthologous genes
  "--orth_genes",
  default = "/CBL_scRNAseq/results/integrated/hgnc_mgi_orth_genes.csv"
)

args <- parser$parse_args()
# for testing
# args <- parser$parse_args(c(
#   "--srat_rds", "/CBL_scRNAseq/results/integrated/20230126/without_future/aldinger_vladoiu_cca.rds",
#   "--out_dir", "/CBL_scRNAseq/results/integrated/20230129"
# ))

if (is.null(args$srat_rds)) {
  stop("Argument for `srat_rds` is missing; please provide an input Seurat object file")
} else {
  srat_rds <- args$srat_rds
}

if (is.null(args$out_dir)) {
  stop("Argument for `out_dir` is missing; please provide an output directory")
} else {
  out_dir <- args$out_dir
  if (!dir.exists(out_dir)) {
    dir.create(out_dir)
  }
}

# ------------------------------------------------------------------------------
# import functions

source("/CBL_scRNAseq/software/utilities/cluster_barplot.R")
source("/CBL_scRNAseq/software/mouse/Vladoiu/add_annotations.R")

# ------------------------------------------------------------------------------
# import data

message(sprintf("Reading in from: %s", srat_rds))
integ_srat <- readRDS(srat_rds)

message(sprintf("Saving output to: %s", out_dir))

# ------------------------------------------------------------------------------
# how many human/mouse cells are in each cluster?

walk2(
  .x = c("stack", "fill"),
  .y = c("num_cells_per_cluster.png", "prop_cells_per_cluster.png"),
  .f = function(x, y, object, split.by, group.by, out_dir) {
    
    # function from `cluster_barplot.R`
    plt <- cluster_barplot(
      object = object,
      split.by = split.by,
      group.by = group.by,
      position = x
    )
    
    # save plot
    ggsave(
      filename = y,
      plot = plt,
      path = out_dir,
      width = 10,
      height = 5,
      units = "in",
      dpi = 600
    )
  },
  object = integ_srat,
  split.by = "species",
  group.by = "seurat_clusters",
  out_dir = out_dir
)

# ------------------------------------------------------------------------------
# what marker genes are expressed in each cluster?

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
m2h_genes <- read_csv(args$orth_genes)
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
plt <- DotPlot(integ_srat, features = genes, assay = "RNA") + 
  scale_x_discrete(limits = rev) + 
  coord_flip()
ggsave(
  "gene_expression_dotplot.png",
  plot = plt,
  path = out_dir,
  width = 12,
  height = 8.5,
  units = "in",
  dpi = 600
)

# ------------------------------------------------------------------------------
# based on Liam's annotations, what human cells are in each cluster?

human_annot <- table(integ_srat$figure_clusters, integ_srat$seurat_clusters) %>% 
  as.data.frame(.)
colnames(human_annot) <- c("annot", "cluster", "freq")

# plot number of annotated human cells in each cluster
plt <- ggplot(human_annot, aes(x = cluster, y = freq, fill = annot)) + 
  geom_col() + 
  labs(y = "number") + 
  theme_light() + 
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) + 
  theme(legend.position = "bottom") + 
  guides(fill = guide_legend(ncol = 7))
ggsave(
  "human_cells_per_cluster.png",
  plot = plt,
  path = out_dir,
  width = 12,
  height = 8,
  units = "in",
  dpi = 600
)

# ------------------------------------------------------------------------------
# based on annotations from Vladoiu paper, what mouse cells are in each cluster?

# change mouse cell types to factor and rename
# function from `add_annotations.R`
integ_srat <- label_vladoiu_cells(integ_srat)

mouse_annot <- table(integ_srat$mouse_cell_type, integ_srat$seurat_clusters) %>% 
  as.data.frame(.)
colnames(mouse_annot) <- c("annot", "cluster", "freq")

# plot number of annotated mouse cells in each cluster
plt <- ggplot(mouse_annot, aes(x = cluster, y = freq, fill = annot)) +
  geom_col() +
  labs(y = "number") +
  theme_light() +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  # use dittoSeq colour palette:
  # https://github.com/dtm2451/dittoSeq/blob/master/R/dittoColors.R
  scale_fill_manual(values = c(
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
  )) + 
  theme(legend.position = "bottom") + 
  guides(fill = guide_legend(ncol = 8))
ggsave(
  "mouse_cells_per_cluster.png",
  plot = plt,
  path = out_dir,
  width = 20,
  height = 10,
  units = "in",
  dpi = 600
)

# ------------------------------------------------------------------------------

print(sessionInfo())
