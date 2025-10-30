# ==============================================================================
# Use `slingshot` to infer pseudotime.
# ==============================================================================

# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(argparse)
library(Seurat)
library(slingshot)
library(tidyverse)
library(patchwork)

# parse command line arguments
parser <- ArgumentParser()

# parser$add_argument(
#   # absolute paths to Seurat RDS files
#   "--srat_rds",
#   default = NULL,
#   required = TRUE
# )
parser$add_argument(
  # output directory
  "--out_dir",
  default = NULL,
  required = TRUE
)

args <- parser$parse_args()

if (is.null(args$out_dir)) {
  stop("Argument for `out_dir` is missing; please provide an output directory")
} else {
  out_dir <- args$out_dir
}

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
# show where human-specific UBCs are clustering

# for harmony, UBCs are clusters 5 (conserved) and 16 (non-homologous/human-specific)
plt <- DimPlot(
  liam_srat,
  reduction = "umap",
  group.by = "integ_clusters",
  cells.highlight = list(
    liam_srat@meta.data[liam_srat$integ_clusters == 16, ] %>% row.names,
    liam_srat@meta.data[liam_srat$integ_clusters == 5, ] %>% row.names
  ),
  # cols.highlight = RColorBrewer::brewer.pal(3, "Set1"),
  sizes.highlight = 0.5
) + 
  scale_colour_manual(
    labels = c("non-UBCs", "5 - homologous", "16 - non-homologous"),
    values = c("grey", RColorBrewer::brewer.pal(3, "Set1")[1:2])
  )

ggsave(
  filename = "non_homologous_ubcs.png",
  plot = plt,
  path = out_dir,
  width = 6,
  height = 4,
  units = "in",
  dpi = 600
)

# ------------------------------------------------------------------------------
# run `slingshot`

# slingshot parameters to loop through
slingshot_args <- bind_rows(
  "pca" = expand_grid(
    cluster_labels = c("new_cell_type", "SCT_snn_res.0.4", "integ_clusters"),
    num_pcs = c(10, 20, 40, Inf)
  ),
  "umap" = expand_grid(
    cluster_labels = c("new_cell_type", "SCT_snn_res.0.4", "integ_clusters"),
    num_pcs = NULL
  ),
  .id = "reduction"
) %>% 
  mutate(start_cluster = case_when(
    cluster_labels == "new_cell_type" ~ "RL-VZ",
    cluster_labels == "SCT_snn_res.0.4" ~ "10",
    cluster_labels == "integ_clusters" ~ "14",
    TRUE ~ NA_character_
  ))

# run slingshot, return list of PseudotimeOrdering objects
all_ptos <- purrr::pmap(
  .l = slingshot_args,
  .f = function(
    srat, 
    reduction, 
    cluster_labels, 
    start_cluster, 
    num_pcs
  ) {
    # run slingshot with error catching
    pto <- tryCatch(
      expr = {
        # custom slingshot helper function from `slingshot.R`
        run_slingshot(
          srat,
          reduction = reduction,
          cluster_labels = cluster_labels,
          start_cluster = start_cluster,
          num_pcs = num_pcs
        )
      },
      error = function(c) {
        # print error message
        message(sprintf(paste0(
          "Error when running slingshot\n",
          c, "\n", # c is the error message that's returned
          "reduction: ", reduction, "\n",
          "cluster_labels: ", cluster_labels, "\n",
          "start_cluster: ", start_cluster, "\n",
          "num_pcs: ", num_pcs, "\n"
        )))
        
        return(NA)
      }
    )
    
    return(pto)
  },
  srat = liam_srat
)

# add names to the pto objects, filtering out the errors
names(all_ptos) <- paste(
  slingshot_args$reduction,
  slingshot_args$num_pcs,
  slingshot_args$cluster_labels,
  slingshot_args$start_cluster,
  sep = "."
)
all_ptos <- all_ptos[!is.na(all_ptos)]


# ------------------------------------------------------------------------------
# save object and make plots

# create directories and save pto object
new_subdirs <- file.path(out_dir, names(all_ptos))
for (subdir in new_subdirs) {
  if (!dir.exists(subdir)) {
    dir.create(subdir)
  }
  message(sprintf("Saving %s", file.path(subdir, "pto.rds")))
  saveRDS(all_ptos[[basename(subdir)]], file = file.path(subdir, "pto.rds"))
}

# colour by cluster
plt_cluster <- purrr::map2(
  .x = all_ptos,
  .y = names(all_ptos),
  .f = function(
    pto,
    subdir
  ) {
    out_subdir <- file.path(out_dir, subdir)
    
    # custom plotting function from `slingshot.R`
    plt <- plot_slingshot(pto, colour_by = "cluster") + 
      labs(title = subdir)
    
    # save as individual file
    message(sprintf("Saving output to %s", out_subdir))
    ggsave(
      "cluster_labels.png",
      plot = plt,
      path = out_subdir,
      width = 6,
      height = 4,
      units = "in",
      dpi = 600
    )
    
    return(plt)
  }
)

# save all to one file
ggsave(
  "all_cluster_labels.png",
  plot = wrap_plots(plt_cluster, ncol = 1),
  path = out_dir,
  width = 6,
  height = 4 * length(plt_cluster),
  units = "in",
  dpi = 600
)

# colour by pseudotime
plt_pseudotime <- purrr::map2(
  .x = all_ptos,
  .y = names(all_ptos),
  .f = function(
    pto,
    subdir
  ) {
    out_subdir <- file.path(out_dir, subdir)
    
    # custom plotting function from `slingshot.R`
    plt <- plot_slingshot(pto, colour_by = "pseudotime")
    
    # save as individual file
    ggsave(
      "pseudotime.png",
      plot = wrap_plots(plt, ncol = 2),
      path = out_subdir,
      width = 10,
      height = 4 * ceiling(length(plt) / 2),
      units = "in",
      dpi = 600
    )
    
    return(plt)
  }
)
