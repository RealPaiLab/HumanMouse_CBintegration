# ==============================================================================
# Integrate Aldinger (human) and Vladoiu (mouse) datasets.
# ==============================================================================

# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# library(biomaRt)
library(argparse)
library(tidyverse)
library(Seurat)

# parse command line arguments
parser <- ArgumentParser()

parser$add_argument(
  # output directory
  "--out_dir",
  default = NULL,
  required = TRUE
)
parser$add_argument(
  # number of cores for t-SNE (and for future parallelization if --future is used)
  "--num_cores",
  default = NULL,
  required = TRUE,
  type = "integer"
)
parser$add_argument(
  # use `future` package for parallelization
  "--future",
  action = "store_true"
)
parser$add_argument(
  # CSV file containing orthologous genes
  "--orth_genes",
  default = "/CBL_scRNAseq/results/integrated/hgnc_mgi_orth_genes.csv"
)
parser$add_argument(
  # name of the file to save the Seurat object
  "--srat_rds",
  default = "aldinger_vladoiu_cca.rds"
)

args <- parser$parse_args()
# args <- parser$parse_args(c("--out_dir", "/CBL_scRNAseq/results/integrated/20230126/", "--num_cores", 24))

if (is.null(args$out_dir)) {
  stop("Argument for `out_dir` is missing; please provide an output directory")
} else {
  out_dir <- args$out_dir
  if (!dir.exists(out_dir)) {
    dir.create(out_dir)
  }
}
num_cores <- as.integer(args$num_cores)
srat_rds <- args$srat_rds

# parallelization with `future` package
# https://satijalab.org/seurat/articles/future_vignette.html
if (args$future) {
  future::plan(strategy = "multisession", workers = num_cores)
  message(sprintf("Using `future` package with %s workers", num_cores))
  # increase RAM available to `future` to 120 GB
  options(future.globals.maxSize = Inf)
}
print(future::plan())

# ------------------------------------------------------------------------------
# import functions

source(file.path("/CBL_scRNAseq/software/utilities/convert_genes.R"))
source(file.path("/CBL_scRNAseq/software/mouse/Vladoiu/add_annotations.R"))

# ------------------------------------------------------------------------------
# import data

aldi_srat <- readRDS("/isilon/CBL_scRNAseq-archived/data/human/Aldinger/seurat.rds")
vlad_srat <- readRDS("/CBL_scRNAseq/results/mouse/Vladoiu/merged_seurat.rds")

# ------------------------------------------------------------------------------
# convert mouse genes to human genes

# get orthologous genes from file
m2h_genes <- read.csv(args$orth_genes)

# remove all duplicated genes
# also remove Pisd (some wonky stuff is happening with that gene)
m2h_genes <- m2h_genes[duplicated(m2h_genes$MGI.symbol) == FALSE
                       & duplicated(m2h_genes$HGNC.symbol) == FALSE
                       & m2h_genes$MGI.symbol != "Pisd", ]

# rename mouse genes to their orthologous human genes
vlad_srat_hgene <- rename_genes(vlad_srat,
                                old_genes = m2h_genes$MGI.symbol,
                                new_genes = m2h_genes$HGNC.symbol)

# ------------------------------------------------------------------------------
# integrate datasets using CCA

obj_list <- list(aldi_srat, vlad_srat_hgene)

# free up some RAM
rm(aldi_srat, vlad_srat, vlad_srat_hgene)
gc()

# re-run sctransform
obj_list <- lapply(
  X = obj_list,
  FUN = SCTransform,
  variable.features.n = 5000,
  vars.to.regress = "CC.Difference",
  return.only.var.genes = TRUE
)

# prep for integration
features <- SelectIntegrationFeatures(obj_list, nfeatures = 5000)
obj_list <- PrepSCTIntegration(obj_list, anchor.features = features)

# find anchors and integrate
anchors <- FindIntegrationAnchors(
  obj_list,
  normalization.method = "SCT",
  anchor.features = features,
  reduction = "cca"
)
integ_srat <- IntegrateData(anchorset = anchors, normalization.method = "SCT")

# normalize "RNA" assay for downstream visualization
message("Normalizing the RNA assay for downstream visualization")
integ_srat <- NormalizeData(integ_srat, assay = "RNA")

message(sprintf(
  "Saving integrated Seurat object as %s. Dimensional reduction not run yet.",
  file.path(out_dir, srat_rds)
))
saveRDS(integ_srat, file = file.path(out_dir, srat_rds))

# free up some RAM
rm(obj_list)
gc()

# ------------------------------------------------------------------------------
# clustering and dimensional reduction

# run PCA, check PCs
message("Running PCA")
integ_srat <- RunPCA(integ_srat, features = VariableFeatures(integ_srat), npcs = 100)
DimHeatmap(integ_srat, dims = 1:15, cells = 500)
ElbowPlot(integ_srat, ndims = 100)

# set ndims for clustering and dimensional reductions (based on elbow plot)
ndims <- 30

# cluster cells
integ_srat <- FindNeighbors(integ_srat, dims = 1:ndims) %>% 
  FindClusters(.)

# run t-SNE
message("Running t-SNE")
integ_srat <- RunTSNE(integ_srat, dims = 1:ndims, num_threads = num_cores)

# run UMAP
message("Running UMAP")
integ_srat <- RunUMAP(integ_srat, dims = 1:ndims)

# ------------------------------------------------------------------------------
# add dataset and species to metadata

integ_srat@meta.data <- integ_srat@meta.data %>% 
  mutate(
    dataset = as.factor(case_when(
      str_detect(orig.ident, "^Vladoiu") ~ "Vladoiu",
      TRUE ~ "Aldinger"
    )),
    species = as.factor(case_when(
      str_detect(orig.ident, "^Vladoiu") ~ "mouse",
      TRUE ~ "human"
    ))
  )

message(sprintf(
  "Saving integrated Seurat object with dimensional reduction results: %s",
  file.path(out_dir, srat_rds)
))
saveRDS(integ_srat, file = file.path(out_dir, srat_rds))

# ------------------------------------------------------------------------------
# save Seurat plots

# make tsne and umap plots, grouped by seurat clusters and by species
for (reduc in c("tsne", "umap")) {
  
  # seurat clusters
  plt1 <- DimPlot(
    integ_srat,
    reduction = reduc,
    group.by = "seurat_clusters",
    label = TRUE,
    repel = TRUE,
    raster = FALSE
  ) + 
    NoLegend()
  
  # species
  plt2 <- DimPlot(
    integ_srat,
    reduction = reduc,
    group.by = "species",
    label = FALSE,
    raster = FALSE
  )
  
  # combine plots and save
  fname <- paste0(reduc, ".png")
  ggsave(
    filename = fname,
    plot = plt1 + plt2,
    path = out_dir,
    width = 16,
    height = 8,
    units = "in",
    dpi = 600
  )
}

# add annotations to Vladoiu dataset (function from `add_annotations.R`)
integ_srat <- label_vladoiu_cells(integ_srat)

# make umap plots grouped by cell types
umap_human_cells <- DimPlot(
  integ_srat,
  cells = rownames(integ_srat@meta.data[integ_srat$species == "human", ]),
  reduction = "umap",
  group.by = "figure_clusters",
  label = TRUE,
  repel = TRUE,
  raster = FALSE
) + 
  ggtitle("human cell types") + 
  theme(legend.position = "bottom") + 
  guides(colour = guide_legend(ncol = 2, override.aes = list(size = 2)))

umap_mouse_cells <- DimPlot(
  integ_srat,
  cells = rownames(integ_srat@meta.data[integ_srat$species == "mouse", ]),
  cols = c(
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
  ),
  reduction = "umap",
  group.by = "mouse_cell_type",
  label = TRUE,
  repel = TRUE,
  raster = FALSE
) + 
  ggtitle("mouse cell types") + 
  theme(legend.position = "bottom") + 
  guides(colour = guide_legend(ncol = 2, override.aes = list(size = 2)))

ggsave(
  filename = "umap_celltypes.png",
  plot = umap_human_cells + umap_mouse_cells,
  path = out_dir,
  width = 18,
  height = 12,
  units = "in",
  dpi = 600
)



# ------------------------------------------------------------------------------

print(sessionInfo())
