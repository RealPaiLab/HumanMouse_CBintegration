# ==============================================================================
# Processes the Vladoiu data and creates a Seurat object using IntegrateData()
# ==============================================================================

# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
# setwd("/CBL_scRNAseq/software/")

library(biomaRt)
library(dplyr)
library(patchwork)
library(Seurat)

# set data and output directories
root_dir <- "CBL_scRNAseq"
data_dir <- file.path("/isilon", root_dir, "data/mouse/Vladoiu")
out_dir <- file.path(paste0("/", root_dir), "results/mouse/Vladoiu")
date_dir <- file.path(out_dir, "20220217")

# ==============================================================================
# function to import data
# ==============================================================================

# load count matrix as Seurat object
load_mtx <- function(
  timepoint, 
  data_dir, 
  mtx = "matrix.mtx.gz", 
  cells = "barcodes.tsv.gz", 
  features = "features.tsv.gz"
) {
  # read in the count matrix
  count_files <- file.path(data_dir, timepoint, c(mtx, cells, features))
  message(sprintf("* Reading in the count matrix for %s", timepoint))
  count_data <- ReadMtx(mtx = count_files[1],
                        cells = count_files[2],
                        features = count_files[3],
                        strip.suffix = TRUE)
  
  # create Seurat object
  message(sprintf("* Creating Seurat object for %s", timepoint))
  srat <- CreateSeuratObject(counts = count_data,
                             project = paste0("Vladoiu-", timepoint),
                             min.cells = 3,
                             min.features = 200)
  
  # rename cell barcodes to something more reasonable
  new_names <- paste0(timepoint, "-", 1:ncol(srat))
  srat <- RenameCells(srat, new.names = new_names)
  
  return(srat)
}

# ==============================================================================
# functions to check QC metrics and filter cells
# ==============================================================================

# plot QC metrics and return plots
plot_qc_metrics <- function(
  srat
) {
  # make violin plot
  qc_plot <- VlnPlot(srat, features = c("nCount_RNA", "nFeature_RNA", "percent_mt"))
  
  # return violin plot
  return(qc_plot)
}

# filter cells
filter_cells <- function (
  srat, 
  deviation = sd, 
  deviation_cutoff = 5
) {
  # calculate percentage of mitochondrial reads
  srat[["percent_mt"]] <- PercentageFeatureSet(srat, pattern = "^mt-")
  
  # calculate the deviation cutoff for total cell count (nCount_RNA) and percent_mt
  count_cut <- deviation(srat[["nCount_RNA"]][, 1]) * deviation_cutoff
  mt_cut <- deviation(srat[["percent_mt"]][, 1]) * deviation_cutoff
  
  srat <- subset(srat, 
                 subset = nCount_RNA < count_cut & percent_mt < mt_cut)
  
  return(srat)
}

# assign cell cycle scores and calculate difference between G2M and S phase
assign_cc_scores <- function(
  srat, 
  s_genes, 
  g2m_genes
) {
  # assign cell cycle score and calculate difference
  srat <- CellCycleScoring(srat, s.features = s_genes, g2m.features = g2m_genes)
  srat$CC.Difference <- srat$S.Score - srat$G2M.Score
  
  return(srat)
}

# ==============================================================================
# functions to integrate list of Seurat objects
# ==============================================================================

integrate_all <- function(
  srat_list, 
  nfeatures = 3000
) {
  features <- SelectIntegrationFeatures(srat_list, nfeatures = nfeatures)
  srat_list <- PrepSCTIntegration(srat_list, anchor.features = features)
  anchors <- FindIntegrationAnchors(srat_list, normalization.method = "SCT", 
                                    anchor.features = features)
  srat_integrated <- IntegrateData(anchors, normalization.method = "SCT")
  return(srat_integrated)
}

# ==============================================================================
# processing the raw count matrices
# ==============================================================================

# get the timepoints
timepoints <- list.dirs(data_dir, full.names = FALSE, recursive = FALSE)

# convert human cell cycling genes to mouse cell cycling genes
if (!exists("s_genes") | !exists("g2m_genes")) {
  # load biomart databases
  source("./utilities/convert_genes.R")
  human <- load_human_bm()
  mouse <- load_mouse_bm()
  
  # convert human genes to mouse genes
  s_genes <- convert_h2m(cc.genes.updated.2019$s.genes, human, mouse)[, 2]
  print(s_genes)
  g2m_genes <- convert_h2m(cc.genes.updated.2019$g2m.genes, human, mouse)[, 2]
  print(g2m_genes)
}

# import count matrices from each timepoint and convert them into a list of Seurat objects
srat_list <- lapply(
  X = timepoints, 
  FUN = function(
    X, 
    data_dir, 
    s_genes,
    g2m_genes
  ) {
    message(sprintf("* Processing the data for timepoint %s", X))
    X <- load_mtx(X, data_dir)
    # need to plot QC and save plots
    X <- filter_cells(X)
    X <- assign_cc_scores(X, s_genes = s_genes, g2m_genes = g2m_genes)
  }, 
  data_dir = data_dir, 
  s_genes = s_genes, 
  g2m_genes = g2m_genes
)

# save list of Seurat objects from the different timepoints
timepoint_list <- file.path(out_dir, "timepoint_seurat_list.rds")
saveRDS(srat_list, file = timepoint_list)

# normalize data with sctransform
srat_list <- lapply(
  X = srat_list, 
  FUN = SCTransform, 
  vars.to.regress = "CC.Difference"
)

# save list of Seurat objects normalized with sctransform
norm_list <- file.path(out_dir, "normalized_seurat_list.rds")
saveRDS(srat_list, file = norm_list)

# integrate data
srat_integ <- integrate_all(srat_list, nfeatures = 3000)

# save integrated object as RDS object
integ_file <- file.path(out_dir, "integrated_seurat.rds")
message(sprintf("Saving the integrated object as %s", integ_file))
saveRDS(srat_integ, file = integ_file)
