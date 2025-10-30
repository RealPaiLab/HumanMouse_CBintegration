# ==============================================================================
# Use `slingshot` and `pseudotimeDE` to infer pseudotime uncertainty.
# ==============================================================================

# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(argparse)
library(Seurat)
library(slingshot)
library(PseudotimeDE)
library(tidyverse)
library(patchwork)
library(parallel)

# parse command line arguments
parser <- ArgumentParser()

parser$add_argument(
  # output directory
  "--out_dir",
  default = NULL,
  required = TRUE
)
parser$add_argument(
  # seed
  "--seed",
  default = 42
)
parser$add_argument(
  # number of subsamples
  "--n_subs",
  default = 1000
)

args <- parser$parse_args()

if (is.null(args$out_dir)) {
  stop("Argument for `out_dir` is missing; please provide an output directory")
} else {
  out_dir <- args$out_dir
}

# set the seed
set.seed(args$seed)

# get the number of subsamples for pseudotime uncertainty
n_subs <- args$n_subs

# import functions
source("/CBL_scRNAseq/software/utilities/slingshot.R")

# ------------------------------------------------------------------------------
# import Seurat objects and copy cluster numbers from the integrated Seurat

liam_srat <- readRDS("/isilon/CBL_scRNAseq/data/human/Aldinger/glutamatergic_dev_Liam.RDS")
integ_srat <- readRDS("/CBL_scRNAseq/results/integrated/20221003/vladoiu_liam_RL_harmony.rds")

# copy the integrated cluster numbers to the human dataset, making sure that the
# cells have the correct cluster labelled
liam_srat$integ_clusters <- integ_srat@meta.data[row.names(liam_srat@meta.data), "seurat_clusters"]

# ------------------------------------------------------------------------------
# run slingshot on all cells, make plots

# custom slingshot function from `slingshot.R`
pto <- run_slingshot(
  srat = liam_srat,
  reduction = "umap",
  cluster_labels = "new_cell_type",
  start_cluster = "RL-VZ"
)

# calculate average pseudotime (`slingshot.R`)
avg_pt <- mean_pseudotime(pto)

# plot average pseudotime
plt_df <- as.data.frame(slingReducedDim(pto))
plt_df$avg_pt <- avg_pt[rownames(plt_df)]
curves <- slingCurves(pto, as.df = TRUE)
plt <- ggplot(plt_df, mapping = aes(x = UMAP_1, y = UMAP_2)) + 
  geom_point(mapping = aes(colour = avg_pt), size = 0.25) + 
  geom_path(data = curves %>% arrange(Order), mapping = aes(group = Lineage)) + 
  theme_classic()
ggsave(
  "mean_pseudotime.png",
  plot = plt,
  path = out_dir,
  width = 5,
  height = 4,
  units = "in",
  dpi = 600
)

# ------------------------------------------------------------------------------
# make subsamples for pseudotimeDE

cell_names <- mclapply(
  X = seq_len(n_subs),
  FUN = function(X, srat) {
    # get the cell names
    cells <- rownames(srat@meta.data)
    
    # sample 80% of cells without replacement
    cells_sample <- sample(x = cells,
                           size = 0.8 * length(cells),
                           replace = FALSE)
    
    return(cells_sample)
  },
  srat = liam_srat
)

# ------------------------------------------------------------------------------
# run slingshot on each of the subsamples for pseudotimeDE



# unfinished: not using pseudotimeDE cuz it doesn't take cell clusters into
# account when calculating differential expression