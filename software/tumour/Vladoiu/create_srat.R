# ==============================================================================
# Load Vladoiu MB scRNA-seq samples from mtx file and create list of Seurat
# objects (one object for each sample).
# ==============================================================================

# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(argparse)
library(tidyverse)
library(Seurat)

# parse command line arguments
parser <- ArgumentParser()

parser$add_argument(
  # input directory with the mtx/barcodes/genes
  "--in_dir",
  default = NULL,
  required = TRUE
)
parser$add_argument(
  # output directory
  "--out_dir",
  default = NULL,
  required = TRUE
)

if (interactive()) {
  # for testing and troubleshooting
  args <- parser$parse_args(c(
    "--in_dir", "/isilon/CBL_scRNAseq-archived/data/tumour/Vladoiu/",
    "--out_dir", "/CBL_scRNAseq/results/tumour/Vladoiu/20230503/"
  ))
} else {
  args <- parser$parse_args()
}

message(sprintf("Saving files to %s", args$out_dir))
if (!dir.exists(args$out_dir)) {
  dir.create(args$out_dir, recursive = TRUE)
}

# ------------------------------------------------------------------------------
# read 10X files (`matrix.mtx`, `barcodes.tsv`, `genes.tsv`)

# directory of each sample
sample_dirs <- list.dirs(args$in_dir, recursive = FALSE) %>% 
  .[str_detect(string = ., pattern = "SHH|G3|G4")]

# make list of Seurat objects (with QC)
srat_list <- purrr::map(
  .x = sample_dirs,
  .f = \(x) {
    message(sprintf("Creating Seurat object for %s", basename(x)))
    srat <- Read10X(file.path(x, "filtered_gene_matrix"), strip.suffix = TRUE) %>% 
      CreateSeuratObject(project = basename(x), min.cells = 3, min.features = 200)
    
    # get % of counts that are from mitochondrial genes
    srat[["percent.mt"]] <- PercentageFeatureSet(srat, pattern = "^MT-")
    
    # set QC cutoffs as per Vladoiu et al. Nature 2019; cells with fewer than
    # 200 genes & genes expressed in fewer than 3 cells were already filtered
    # out in `CreateSeuratObject`
    count_max <- median(srat$nCount_RNA) + 4*sd(srat$nCount_RNA)
    mito_max <- median(srat$percent.mt) + 4*sd(srat$percent.mt)
    
    # filter cells
    srat <- subset(srat, subset = nCount_RNA < count_max & percent.mt < mito_max)
    
    return(srat)
  }
)

# ------------------------------------------------------------------------------
# normalize with sctransform (each sample separately)

# first, run cell cycle scoring
s_genes <- cc.genes$s.genes
g2m_genes <- cc.genes$g2m.genes

srat_list <- purrr::map(
  .x <- srat_list,
  .f = \(x) {
    # cell cycle scoring
    x <- CellCycleScoring(x, s.features = s_genes, g2m.features = g2m_genes)
    x$CC.Difference <- x$S.Score - x$G2M.Score
    
    # regress out the cell cycle scoring
    message(sprintf("Running `SCTransform` on %s", unique(x$orig.ident)))
    x <- SCTransform(
      object = x,
      vars.to.regress = "CC.Difference",
      variable.features.n = NULL
    )
  }
)

# ------------------------------------------------------------------------------
# merge samples together

srat_merged <- merge(
  srat_list[[1]],
  srat_list[-1],
  add.cell.ids = basename(sample_dirs),
  project = "vladoiu_mb"
)

# set variable features (https://github.com/satijalab/seurat/issues/5135)
# DO NOT USE `FindVariableFeatures()` cuz it cannot be used with SCTransform
VariableFeatures(srat_merged[["SCT"]]) <- rownames(srat_merged[["SCT"]]@scale.data)

# ------------------------------------------------------------------------------
# run dimensionality reduction

message("Running PCA")
srat <- RunPCA(srat_merged, npcs = 100)

ndims <- 50

message("Running UMAP")
srat <- RunUMAP(srat, dims = 1:ndims)

srat <- FindNeighbors(srat, dims = 1:ndims) %>% 
  FindClusters()

# ------------------------------------------------------------------------------
# make plots

# elbow plot
ggsave(
  filename = "pca_elbow.png",
  plot = ElbowPlot(srat, ndims = 100),
  path = args$out_dir,
  width = 8,
  height = 6,
  units = "in",
  dpi = 600
)

# PCA
plt <- DimPlot(srat, group.by = "orig.ident", reduction = "pca") + 
  DimPlot(srat, group.by = "seurat_clusters", reduction = "pca", label = TRUE)
ggsave(
  filename = "pca.png",
  plot = plt,
  path = args$out_dir,
  width = 15,
  height = 6,
  units = "in",
  dpi = 600
)

# UMAP
plt <- DimPlot(srat, group.by = "orig.ident") + 
  DimPlot(srat, group.by = "seurat_clusters", label = TRUE)
ggsave(
  filename = "umap.png",
  plot = plt,
  path = args$out_dir,
  width = 15,
  height = 6,
  units = "in",
  dpi = 600
)

# ------------------------------------------------------------------------------
# save to file

saveRDS(
  object = srat,
  file = file.path(args$out_dir, "vladoiu_mb_merge.rds")
)

# ------------------------------------------------------------------------------

message("\nSESSION INFO\n")
print(sessionInfo())

