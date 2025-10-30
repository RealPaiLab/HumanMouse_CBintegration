# ==============================================================================
# Calculate and plot the difference in species centroids grouped by cell types.
# ==============================================================================

# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(argparse)
library(tidyverse)
library(Seurat)

# parse command line arguments
parser <- ArgumentParser()
# absolute paths to Seurat RDS files (can pass multiple files)
parser$add_argument(
  "--srat_rds",
  action = "extend",
  nargs = "+",
  default = NULL,
  required = TRUE
)
# integration method for each Seurat RDS file listed (e.g., CCA, RPCA, harmony, etc.)
parser$add_argument(
  "--integ_method",
  action = "extend",
  nargs = "+",
  default = NULL,
  required = TRUE
)
parser$add_argument(
  # output directory
  "--out_dir",
  default = NULL,
  required = TRUE
)

args <- parser$parse_args()

# srat_rds and integ_method are required arguments
srat_rds <- args$srat_rds
integ_method <- args$integ_method

# make sure number of files and method labels are of same length
if (length(srat_rds) != length(integ_method)) {
  stop("srat_rds and integ_method must have same number of elements")
}

if (is.null(args$out_dir)) {
  stop("Argument for `out_dir` is missing; please provide an output directory")
} else {
  out_dir <- args$out_dir
  if (!dir.exists(out_dir)) {
    dir.create(out_dir)
  }
}

message(sprintf("Saving results to: %s", out_dir))

# import functions
source(file.path("/CBL_scRNAseq/software/mouse/Vladoiu/add_annotations.R"))
source(file.path("/CBL_scRNAseq/software/utilities/cell_labelling.R"))
source(file.path("/CBL_scRNAseq/software/utilities/score_integration.R"))

# ------------------------------------------------------------------------------
# calculate distance between centroids for all items

# import Seurat objects and store as list
all_srat <- purrr::map(
  .x = srat_rds,
  .f = readRDS
)
names(all_srat) <- integ_method

# get distance between human and mouse centroids
all_centroid_dist <- purrr::map2_dfr(
  .x = all_srat,
  .y = names(all_srat),
  .f = function(x, y) {
    # add mouse cell types (`add_annotation.R`)
    x <- label_vladoiu_cells(x) %>% 
      # pool cell types (`cell_labelling.R`)
      pool_cell_types(col_name = "common_cell_type")
    
    # get embeddings and metadata
    if ("harmony" %in% names(x@reductions)) {
      embeddings <- Embeddings(x, reduction = "harmony")
    } else {
      embeddings <- Embeddings(x, reduction = "pca")
    }
    metadata <- x@meta.data
    
    # calculate distance between human and mouse centroids for each cell type
    # (`score_integration.R`)
    all_dist <- delta_centroids(
      embeddings = embeddings,
      metadata = metadata,
      by = "common_cell_type",
      filter_out = c("other/missing"),
      integ_method = y
    )

    return(all_dist)
  }
)

# ------------------------------------------------------------------------------
# plot delta centroids

plt <- ggplot(all_centroid_dist, aes(x = common_cell_type, y = delta)) + 
  geom_point(aes(colour = integ_method, shape = integ_method), size = 3) + 
  labs(y = "distance between mouse\nand human cluster centroids") + 
  scale_y_continuous(limits = c(0, NA)) + # `NA` uses the current min/max
  scale_color_brewer(palette = "Set1") + 
  theme_light() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

ggsave(
  filename = "delta_centroids.png",
  plot = plt,
  path = out_dir,
  width = 5,
  height = 4,
  units = "in",
  dpi = 600
)

# ------------------------------------------------------------------------------
# plot distribution of the embeddings

all_embeds <- map2_dfr(
  .x = all_srat,
  .y = names(all_srat),
  .f = function(x, y) {
    # get embeddings
    if ("harmony" %in% names(x@reductions)) {
      embeddings <- Embeddings(x, reduction = "harmony")
    } else {
      embeddings <- Embeddings(x, reduction = "pca")
    }
    
    # coerce to data frame
    embeddings <- as.data.frame(embeddings) %>%
      # add rownames as a column
      rownames_to_column(var = "cell_id") %>%
      # pivot longer
      pivot_longer(cols = !cell_id,
                   names_to = "pc",
                   values_to = "embedding") %>%
      # add integration method as column, remove prefix in the PC column
      mutate(integ_method = y,
             pc = str_replace(pc, "^[:alpha:]+", "PC"))
    
    embeddings$pc <- factor(embeddings$pc, levels = str_sort(unique(embeddings$pc), numeric = TRUE))
    
    return(embeddings)
  }
)

plt <- ggplot(all_embeds, aes(x = integ_method, y = embedding, colour = integ_method)) + 
  geom_boxplot(outlier.size = 1) + 
  facet_wrap(vars(pc), ncol = 10, scales = "free") + 
  scale_colour_brewer(palette = "Set1")

ggsave(
  filename = "principal_components_boxplots.png",
  plot = plt,
  path = out_dir,
  width = 20,
  height = 20,
  units = "in",
  dpi = 600
)

plt <- ggplot(all_embeds, aes(x = embedding, colour = integ_method)) + 
  geom_density(aes(fill = integ_method), alpha = 0.1) + 
  facet_wrap(vars(pc), ncol = 10, scales = "free") + 
  scale_colour_brewer(palette = "Set1") + 
  scale_fill_brewer(palette = "Set1")

ggsave(
  filename = "principal_components_density.png",
  plot = plt,
  path = out_dir,
  width = 20,
  height = 20,
  units = "in",
  dpi = 600
)

# ------------------------------------------------------------------------------

print(sessionInfo())
