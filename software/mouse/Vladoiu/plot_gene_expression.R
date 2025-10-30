# ==============================================================================
# Take the Vladoiu data and plot the expression of various marker genes
# ==============================================================================

# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
# setwd("/CBL_scRNAseq/software/")

library(tidyverse)
library(Seurat)

# set data and output directories
root_dir <- "CBL_scRNAseq"
data_dir <- file.path("/isilon", root_dir, "data/mouse/Vladoiu")
out_dir <- file.path("", root_dir, "results/mouse/Vladoiu")
date_dir <- file.path(out_dir, format(Sys.Date(), "%Y%m%d"))

# create output directory if it doesn't exist
if (!dir.exists(date_dir)) {
  dir.create(date_dir)
}

# load Seurat object from RDS file
srat <- readRDS(file.path(out_dir, "merged_seurat.rds"))

# ==============================================================================
# functions to make and save violin/feature plots
# ==============================================================================

# make violin plot and save it
save_vln <- function(
  srat,
  gene, 
  group.by = NULL, 
  slot = NULL, 
  show_plot = FALSE, 
  path = NULL, 
  width = 11.5, 
  height = 8.5
) {
  message(sprintf("* Making violin plot of %s", gene))
  
  # make violin plot
  vln <- VlnPlot(srat, features = gene, group.by = group.by, slot = slot %||% "data")
  
  # show plots
  if (show_plot) {
    print(vln)
  }
  
  # create filename
  filename <- paste0(gene, "_vln.pdf")
  
  # save violin plot
  ggsave(filename = filename, plot = vln, path = path,
         width = width, height = height, units = "in")
}

# make feature plots and save it
save_feature <- function (
  srat, 
  gene, 
  order = FALSE, 
  slot = NULL, 
  min.cutoff = NA, 
  max.cutoff = NA, 
  show_plot = FALSE, 
  path = NULL, 
  width = 5, 
  height = 4
) {
  message(sprintf("* Making feature plot of %s", gene))
  
  # make feature plots for t-SNE and UMAP
  for (dimred in c("tsne", "umap")) {
    plt <- FeaturePlot(
      srat, 
      features = gene, 
      reduction = dimred, 
      order = order, 
      slot = slot %||% "data", 
      min.cutoff = min.cutoff, 
      max.cutoff = max.cutoff
    )
    
    # show plots
    if (show_plot) {
      print(plt)
    }
    
    # create filename
    filename <- paste0(gene, "_", dimred, ".pdf")
    
    # save feature plot
    ggsave(filename = filename, plot = plt, path = path, 
           width = width, height = height, units = "in")
  }
}

# ==============================================================================
# save plots
# ==============================================================================

# major gene markers from Carter et al. 2018 and Vladoiu et al 2019
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

# set path to save plots
plot_path <- file.path(date_dir, "gene_expression_plots")
if (!dir.exists(plot_path)) {
  dir.create(plot_path)
}

# loop through genes to generate plots
for (gene in genes) {
  save_vln(srat, gene = gene, path = plot_path)
  save_feature(srat, gene = gene, order = TRUE, 
               min.cutoff = "q10", max.cutoff = "q90", 
               path = plot_path)
}

# save a DotPlot of the same genes
dot <- DotPlot(srat, features = genes) + scale_x_discrete(limits = rev) + coord_flip()
ggsave("gene_expression_dotplot.pdf", plot = dot, path = plot_path, 
       width = 11.5, height = 8.5, units = "in")
