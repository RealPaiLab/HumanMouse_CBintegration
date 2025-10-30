# ==============================================================================
# Processes the Vladoiu data and creates a Seurat object by merging the timepoints
# ==============================================================================

# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
# setwd("/CBL_scRNAseq/software/")

library(tidyverse)
library(Seurat)

# set data and output directories
root_dir <- "CBL_scRNAseq"
data_dir <- file.path("/isilon", root_dir, "data/mouse/Vladoiu")
out_dir <- file.path("", root_dir, "results/mouse/Vladoiu")

# ==============================================================================
# import functions 
# ==============================================================================

source(file.path("", root_dir, "software/mouse/Vladoiu/load_vladoiu.R"))

# ==============================================================================
# load and merge the count matrices from the 9 timepoints
# ==============================================================================

# get the timepoints
timepoints <- list.dirs(data_dir, full.names = FALSE, recursive = FALSE)

# load count matrices into a list of Seurat objects
# filter cells with high counts or high mitochondrial content
srat_list <- lapply(
  X = timepoints, 
  FUN = function(X) {
    X <- load_mtx(X, data_dir)
    X <- filter_cells(X, deviation = sd, deviation_cutoff = 5)
  }
)

# merge counts into one Seurat object
srat_merge <- merge(
  x = srat_list[[1]], 
  y = srat_list[-1], 
  project = "Vladoiu"
)

# ==============================================================================
# process and normalize the merged data
# ==============================================================================

# convert human cell cycling genes to mouse cell cycling genes
if (!exists("s_genes") | !exists("g2m_genes")) {
  source(file.path("", root_dir, "software/utilities/convert_genes.R"))
  
  # load orthologous genes from file
  m2h_genes <- read.csv(file.path("", root_dir, "results/integrated/hgnc_mgi_orth_genes.csv"))
  
  # convert genes
  s_genes <- get_orth_genes(m2h_genes$HGNC.symbol, m2h_genes$MGI.symbol, cc.genes.updated.2019$s.genes)
  g2m_genes <- get_orth_genes(m2h_genes$HGNC.symbol, m2h_genes$MGI.symbol, cc.genes.updated.2019$g2m.genes)
}

# assign cell cycle score to each cell
srat_merge <- CellCycleScoring(srat_merge, s.features = s_genes, g2m.features = g2m_genes)
srat_merge$CC.Difference <- srat_merge$S.Score - srat_merge$G2M.Score

# normalize with sctransform
srat_merge <- SCTransform(srat_merge, vars.to.regress = "CC.Difference")

# add actual annotations from Vladoiu paper
source("./add_annotations.R")
srat_merge <- add_annotations(srat_merge)

# save merged Seurat as RDS
merge_file <- file.path(out_dir, "merged_seurat.rds")
message(sprintf("Saving the merged object as %s", merge_file))
saveRDS(srat_merge, file = merge_file)
