# ==============================================================================
# Generate QC plots for Taylor lab MB tumour scRNA-seq samples and combine into
# one Seurat object.
# ==============================================================================

library(argparse)
library(tidyverse)
library(patchwork)
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
  arg_list <- parser$parse_args(c(
    "--in_dir", "/.mounts/labs/pailab/private/projects/MB_multiome",
    "--out_dir", "/.mounts/labs/pailab/private/icheong/CBL_scRNAseq/results/tumour/MB_multiome/20250204"
  ))
} else {
  arg_list <- parser$parse_args()
}

message(sprintf("Saving files to %s", getwd()))
if (!dir.exists(arg_list$out_dir)) {
  dir.create(arg_list$out_dir, recursive = TRUE)
}

# ------------------------------------------------------------------------------
# functions

plot_qc_metrics <- function(
  srat,
  filename,
  path
) {
  if (!dir.exists(path)) {
    dir.create(path, recursive = TRUE)
  }
  # plot number of genes, number of counts, percent mitochondrial genes
  .plt <- map(
    .x = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
    .f = \(x) {
      stats <- sprintf("range: %s-%s", round(min(srat[[x]]), 2), round(max(srat[[x]]), 2))
      .plt <- VlnPlot(srat, features = x, group.by = "sample_id") + 
        NoLegend() + 
        labs(caption = stats) + 
        theme(
          axis.text.x = element_text(angle = 0, hjust = 0.5)
        )
      return(.plt)
    }
  )
  ggsave(
    filename = filename,
    plot = wrap_plots(.plt),
    path = path,
    width = 10,
    height = 6,
    units = "in",
    dpi = 600
  )
}

plot_umaps <- function(
  srat,
  filename,
  path
) {
  if (!dir.exists(path)) {
    dir.create(path, recursive = TRUE)
  }
  # plot previously generated UMAPs
  .plt <- map(
    .x = c("umap.rna", "umap.atac", "umap"),
    .f = \(reduc) {
      .plt <- DimPlot(
        srat,
        reduction = reduc,
        group.by = "seurat_clusters",
        label = FALSE
      ) + 
        labs(title = reduc)
    }
  )
  ggsave(
    filename = filename,
    plot = wrap_plots(.plt, guides = "collect"),
    path = path,
    width = 12,
    height = 4,
    units = "in",
    dpi = 600
  )
}

# ------------------------------------------------------------------------------
# load sample metadata

# sample_md <- read_tsv(
#   file = file.path(arg_list$in_dir, "20241212_multiome_metadata.tsv")
# ) %>%
#   mutate(
#     subtype = str_extract(
#       string = methyl_dx,
#       pattern = "MB, \\b([:alnum:]+)\\b",
#       group = 1
#     )
#   )

# JSON file from Gabrielle Persad (Stein lab)
qc_filters <- jsonlite::fromJSON(file.path(arg_list$in_dir, "MB_qc_filters.json")) %>%
  pluck(1)

# ------------------------------------------------------------------------------
# load each dataset and plot QC metrics

# for each sample:
# - make QC plots for the sample
# - extract QC stats so they can all be plotted on one page
qc_stats <- map(
  .x = names(qc_filters),
  .f = \(sample_id) {
    message(sprintf("***Reading in and generating plots for %s***", sample_id))
    srat <- readRDS(file.path(arg_list$in_dir, "input", paste0(sample_id, ".rds")))

    # QC plots
    plot_qc_metrics(
      srat = srat,
      filename = paste0(sample_id, "_qc_metrics.png"),
      path = file.path(arg_list$out_dir, "qc_plots")
    )

    # plot previously generated UMAPs
    plot_umaps(
      srat = srat,
      filename = paste0(sample_id, "_orig_umaps.png"),
      path = file.path(arg_list$out_dir, "umaps")
    )

    # get QC stats for the sample
    qc_df <- dplyr::select(
      srat[[]],
      sample_id,
      nFeature_RNA,
      nCount_RNA,
      percent.mt
    ) %>%
      remove_rownames()

    # add subtype
    qc_df$subtype = qc_filters[[sample_id]]$subtype

    return(qc_df)
  }
) %>%
  do.call(what = rbind, args = .) %>%
  mutate(subtype = fct_relevel(subtype, "WNT", "SHH", "G3", "G4/G3", "G4"))

# plot QC metrics for all samples on one plot
qc_plt <- map(
  .x = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  .f = \(metric) {
    .plt <- ggplot(
      qc_stats,
      aes(x = sample_id, y = .data[[metric]], fill = subtype)
    ) + 
      geom_violin(linewidth = 0.25) + 
      facet_grid(cols = vars(subtype), scales = "free", space = "free") + 
      labs(title = metric) + 
      theme_light() + 
      theme(
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        legend.position = "none"
      )

    return(.plt)
  }
)

# plot number of cell in each sample
ncells_plt <- ggplot(qc_stats, aes(sample_id, fill = subtype)) + 
  geom_bar(width = 0.75) + 
  facet_grid(cols = vars(subtype), scales = "free", space = "free") + 
  labs(title = "Number of cells") + 
  theme_light() + 
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1),
    legend.position = "none"
  )

.plt <- append(qc_plt, list(ncells_plt)) %>%
  wrap_plots(ncol = 1)

ggsave(
  filename = "all_samples.png",
  plot = .plt,
  path = file.path(arg_list$out_dir, "qc_plots"),
  width = 24,
  height = 12,
  units = "in",
  dpi = 600
)

# ------------------------------------------------------------------------------

message("\n***SESSION INFO***\n")
print(sessionInfo())

