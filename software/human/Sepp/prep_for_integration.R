# ==============================================================================
# Prepare Sepp 2024 dataset for integration. See QC in 2024-03-26 notes.
# ==============================================================================

# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(argparse)
library(tidyverse)
library(Seurat)
library(SingleCellExperiment)
library(biomaRt)

# parse command line arguments
parser <- ArgumentParser()

parser$add_argument(
  # subset cells for testing
  "--gene_table",
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
  arg_list <- parser$parse_args(c(
    "--gene_table", "/CBL_scRNAseq/results/human/Sepp/20240329/ensg_hgnc_ids.csv",
    "--out_dir", "/CBL_scRNAseq/results/human/Sepp/20240329/"
  ))
} else {
  arg_list <- parser$parse_args()
}

message(sprintf("***Saving files to %s***", arg_list$out_dir))
if (!dir.exists(arg_list$out_dir)) {
  dir.create(arg_list$out_dir, recursive = TRUE)
}

# load functions
source("/CBL_scRNAseq/software/utilities/convert_genes.R")

# ------------------------------------------------------------------------------
# load SCE object and convert to Seurat object

sce <- readRDS("/isilon/CBL_scRNAseq-archived/data/src/neurodev-genomics/scRNAseq/Sepp_2024/hum_sce_final.rds")

srat <- as.Seurat(
  sce,
  counts = "umi",
  data = NULL,
  project = "sepp_human"
) %>% 
  RenameAssays(originalexp = "RNA")

rm(sce) # clean up memory

# ------------------------------------------------------------------------------
# filter out cells collected from postnatal samples

# get prenatal samples
prenatal <- unique(srat$Stage) %>% 
  str_subset(pattern = "wpc$")

# subset cells
num_cells_before <- ncol(srat)
srat <- subset(
  srat,
  subset = Stage %in% prenatal
)
num_cells_after <- ncol(srat)

# print number of cells that were filtered
message(sprintf(
  "***Removed %s postnatal cells; %s prenatal cells remaining***",
  num_cells_before - num_cells_after,
  num_cells_after
))

# ------------------------------------------------------------------------------
# convert Ensembl gene IDs to HGNC symbol

# read in table mapping Ensembl to HGNC genes
gene_table <- read_csv(
  file = arg_list$gene_table
)

# get Ensembl genes that don't have an HGNC symbol or Ensembl genes that have
# multiple HGNC symbols
genes_to_remove <- gene_table$ensembl_gene_id[
  is.na(gene_table$hgnc_symbol) | 
    duplicated(gene_table$hgnc_symbol) | 
    duplicated(gene_table$hgnc_symbol, fromLast = TRUE)
]

# get Ensembl genes that couldn't be found (due to different Ensembl versions)
genes_to_remove <- c(
  genes_to_remove,
  rownames(srat)[!(rownames(srat) %in% gene_table$ensembl_gene_id)]
) %>% 
  unique()

# get genes to keep
genes_to_keep <- rownames(srat)[!(rownames(srat) %in% genes_to_remove)]

# subset genes
srat <- subset(srat, features = genes_to_keep)

# print number of genes that were filtered
message(sprintf(
  "***Removed %s genes; %s genes remaining***",
  length(genes_to_remove),
  length(genes_to_keep)
))

# rename the genes
filtered_gene_table <- gene_table[gene_table$ensembl_gene_id %in% genes_to_keep, ]
srat <- rename_genes(
  object = srat,
  old_genes = filtered_gene_table$ensembl_gene_id,
  new_genes = filtered_gene_table$hgnc_symbol
)

# ------------------------------------------------------------------------------
# remove duplicated metadata columns and set levels for sample age

srat@meta.data <- srat@meta.data %>% 
  dplyr::select(
    -c(nCount_originalexp, nFeature_originalexp)
  ) %>% 
  relocate(
    c(nCount_RNA, nFeature_RNA),
    .after = orig.ident
  )

srat$Stage <- fct_relevel(
  srat$Stage,
  levels = str_sort(unique(srat$Stage), numeric = TRUE)
)

# ------------------------------------------------------------------------------
# calculate percentage of mitochondrial RNA

srat$percent.mt <- PercentageFeatureSet(srat, pattern = "^MT-")
srat@meta.data <- relocate(srat@meta.data, percent.mt, .after = "nFeature_RNA")

# ------------------------------------------------------------------------------
# normalize and regress out cell cycle

# cell cycle scoring
srat <- CellCycleScoring(
  srat,
  s.features = cc.genes.updated.2019$s.genes,
  g2m.features = cc.genes.updated.2019$g2m.genes
)
srat$CC.Difference <- srat$S.Score - srat$G2M.Score

# sctransform
srat <- SCTransform(
  srat,
  variable.features.n = 5000,
  vars.to.regress = "CC.Difference"
)

# ------------------------------------------------------------------------------
# run PCA, UMAP, clustering

srat <- RunPCA(srat, assay = "SCT", npcs = 100)
srat <- RunUMAP(srat, assay = "SCT", dims = 1:25)
srat <- FindNeighbors(srat, dims = 1:25)
srat <- FindClusters(srat, resolution = 0.8)

# ------------------------------------------------------------------------------
# save Seurat object

srat_out_path <- file.path(arg_list$out_dir, "sepp2024_sct.rds")
message(sprintf("\n***Saving processed Seurat object to %s***\n", srat_out_path))
saveRDS(object = srat, file = srat_out_path)

# ------------------------------------------------------------------------------

message("\n***SESSION INFO***\n")
print(sessionInfo())
