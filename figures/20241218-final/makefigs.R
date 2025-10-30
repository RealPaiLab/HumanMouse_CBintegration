# ==============================================================================
# These figures were generated for the OICR PI seminar on 2024-09-24. This
# script should be run in the `scrnaseq_env` conda environment.
#
# The EnrichmentMap of the non-orthologous genes (`nonorth_human.svg`) was
# copied from the integrated results folder from 2024-10-18.
# ==============================================================================

library(tidyverse)
library(Seurat)
library(patchwork)
library(ComplexHeatmap)

# import functions
source("/.mounts/labs/pailab/private/icheong/CBL_scRNAseq/software/utilities/cell_labelling.R")
source("/.mounts/labs/pailab/private/icheong/CBL_scRNAseq/software/utilities/plotting.R")
source("/.mounts/labs/pailab/private/icheong/CBL_scRNAseq/thesis/utils.R")

# colour palettes (function from `utils.R`)
my_pals <- get_custom_pals()

# >>> load Seurat objects >>>

# `get_srat_paths` and `load_srat` from `utils.R`
srat_qs <- get_srat_paths()
srat <- load_srat(srat_qs[c("ubc", "mb")])

# <<< load Seurat objects <<<

# ------------------------------------------------------------------------------
# full cerebellum UMAP; code taken from `integ_full.qmd` in thesis figures
# directory

# load Seurat object for integrated dataset
full_srat <- load_srat(srat_qs["full"]) %>%
  pluck("full")

# plot UMAP by species
.plt <- DimPlot(
  full_srat,
  reduction = "umap",
  cells.highlight = WhichCells(full_srat, expression = species == "human"),
  sizes.highlight = 0.01,
  label = FALSE,
  raster = FALSE
) + 
  scale_colour_manual(
    labels = c("mouse", "human"),
    values = my_pals$species
  ) + 
  theme_classic2() + 
  theme(
    legend.position = "inside",
    legend.position.inside = c(0.025, 0.25),
    legend.justification = c(0, 0.5)
  )
ggsave(
  filename = "full_umap_species.png",
  plot = .plt,
  width = 5,
  height = 5,
  units = "in",
  dpi = 600
)

# plot UMAP by cell type
.plt <- DimPlot(
  full_srat,
  reduction = "umap",
  group.by = "common_cell_name",
  label = TRUE,
  label.size = 3,
  repel = TRUE,
  raster = FALSE
) + 
  labs(title = NULL) + 
  theme_classic2() + 
  theme(legend.position = "none")
ggsave(
  filename = "full_umap_cell_type.png",
  plot = .plt,
  width = 5,
  height = 5,
  units = "in",
  dpi = 600
)

# ------------------------------------------------------------------------------
# pathway enrichment analysis for human-enriched UBCs

# load pathway enrichment results for UBC_1 and UBC_2
ubc_subclusters <- paste0("UBC_", c(1:2))
ubc_pathways <- map(
  .x = ubc_subclusters,
  .f = \(clust) {
    ubc_pathways <- read_csv(file.path(
      "/.mounts/labs/pailab/private/icheong/CBL_scRNAseq/results/integrated/20241029",
      clust,
      "enr_pathways.csv"
    ))

    return(ubc_pathways)
  }
) %>%
  setNames(nm = ubc_subclusters)

walk(
  .x = ubc_subclusters,
  .f = \(clust) {
    # subset pathways for the subcluster
    pway <- ubc_pathways[[clust]]

    # number of significant pathways
    n_sig <- nrow(filter(pway, p_value < 0.05))

    # total number of enriched pathways
    caption <- paste0("# enriched pathways/total = ", n_sig, "/", nrow(pway))

    # how many pathways to show?
    n_plot <- 8

    # get pathways to plot
    pway <- slice_min(
      pway,
      order_by = p_value,
      n = n_plot
    )

    # set font size and width of terms
    term_font_size <- rel(1)
    lab_width <- 25

    # plot enriched pathways
    .plt <- ggplot(
      pway,
      aes(
        x = -log10(p_value),
        y = fct_reorder(.f = term_name, .x = p_value, .desc = TRUE)
      )
    ) + 
      geom_col(fill = my_pals$ubc_integ_clust[[clust]]) + 
      labs(
        caption = caption,
        x = "-log10(FDR)",
        y = "pathway term"
      ) + 
      scale_x_continuous(expand = expansion(mult = c(0, 0.05))) + 
      scale_y_discrete(label = scales::label_wrap(lab_width)) + 
      theme_classic2() + 
      theme(
        axis.text.y = element_text(size = term_font_size),
        plot.title = element_text(vjust = 1)
      )
    ggsave(
      filename = paste0(clust, "_pathways.png"),
      plot = .plt,
      width = 4,
      height = 4.5,
      units = "in",
      dpi = 600
    )
  }
)

# ------------------------------------------------------------------------------
# regulon expression UMAP

# load integrated human cells
human_rl <- readRDS("/.mounts/labs/pailab/private/icheong/CBL_scRNAseq/results/pyscenic/20240704/integ_human_srat.rds")

# load regulon expression
regulon_expr <- read.csv("/.mounts/labs/pailab/private/icheong/CBL_scRNAseq/results/pyscenic/20240926/aldinger_sepp_RL_regulon_expr.csv", row.names = 1) %>%
  select(matches("OTX2|SOX4|SOX11")) %>%
  rename_with(.fn = \(x) {
    str_replace(string = x, pattern = "Regulon\\.(.*)\\.", replacement = "\\1_reg")
  })

# add regulon expression to metadata
human_rl@meta.data <- merge(
  x = human_rl[[]],
  y = regulon_expr,
  by = "row.names",
  all.x = TRUE
) %>%
  column_to_rownames("Row.names")

# plot regulon expression
.plt <- FeaturePlot(
  human_rl,
  features = paste0(c("OTX2", "SOX4", "SOX11"), "_reg"),
  order = TRUE,
  combine = FALSE
) %>%
  setNames(nm = c("OTX2", "SOX4", "SOX11"))

walk(
  .x = names(.plt),
  .f = \(reg) {
    ggsave(
      filename = paste0(reg, "_regulon.png"),
      plot = .plt[[reg]] + theme_classic2(),
      width = 5,
      height = 5,
      units = "in",
      dpi = 600
    )
  }
)

# ------------------------------------------------------------------------------
# Group 4 MB volcano plot (copied from mb_similarity.qmd)

# load list of TFs
tfs <- read.csv(
  "/.mounts/labs/pailab/src/transcription-factors/human-tfs-lambert/full_database_v1.01.csv",
  row.names = 1,
  check.names = FALSE
) %>%
  # get only confirmed TFs
  dplyr::filter(`Is TF?` == "Yes") %>%
  pull(`HGNC symbol`)

# load differential gene expression results for MB subtypes
mb_diff_expr <- read_csv("/.mounts/labs/pailab/private/icheong/CBL_scRNAseq/results/tumour/Vladoiu/20240917/mb_subtype_diff_exp.csv")

# custom volcano plot function (slightly different from the one in `plotting.R`)
make_volcano_plot <- function(
  data,
  log_fc,
  log_pval,
  direction,
  gene_labels
) {
  .plt <- ggplot(
    data = data,
    mapping = aes(
      x = {{ log_fc }},
      y = {{ log_pval }}
    )
  ) + 
    geom_point(
      aes(colour = {{ direction }}),
      alpha = 0.5,
      shape = "circle",
      size = 2,
      stroke = 0
    ) + 
    ggrepel::geom_text_repel(
      aes(label = {{ gene_labels }}, colour = {{ direction }}),
      size = 3,
      min.segment.length = 0.1,
      max.overlaps = Inf,
      show.legend = FALSE,
      seed = 42
    ) + 
    geom_vline(
      xintercept = 0,
      colour = "grey",
      linewidth = 0.25,
      linetype = "dashed"
    ) + 
    geom_hline(
      yintercept = -log10(0.05),
      colour = "grey",
      linewidth = 0.25,
      linetype = "dashed"
    )

  return(.plt)
}

# set genes with `p_val_adj` of 0 to 1e-300 (the minimum non-zero `p_val_adj` is
# around 3e-294 right now)
min_adj_p <- filter(mb_diff_expr, p_val_adj != 0) %>% pull(p_val_adj) %>% min()

# prep for plotting
g4_diff_expr_genes <- map(
  .x = "G4",
  .f = \(subtype) {
    # genes to label
    label_genes <- c(
      "EOMES",
      "OTX2",
      "SOX4",
      "SOX11",
      "FOXP2",
      "CBFA2T2",
      "TBR1",
      "NEUROD1",
      "FOS",
      "JUN"
    )

    de <- filter(mb_diff_expr, cluster == subtype) %>%
      mutate(
        direction = case_when(
          avg_log2FC < 0 & p_val_adj < 0.05 ~ "down",
          avg_log2FC > 0 & p_val_adj < 0.05 ~ "up",
          .default = "n.s."
        ),
        p_val_adj = case_when(
          p_val_adj < min_adj_p ~ 1e-300,
          .default = p_val_adj
        ),
        tf_label = case_when(
          gene %in% label_genes ~ gene,
          .default = ""
        )
      ) %>%
      arrange(tf_label)

    return(de)
  }
)

# volcano plots of TFs only (all subtypes)
walk(
  .x = "G4",
  .f = \(subtype) {
    de <- g4_diff_expr_genes[[1]] %>%
      filter(gene %in% tfs)

    # label number of TFs
    num_up <- nrow(filter(de, direction == "up"))
    num_down <- nrow(filter(de, direction == "down"))
    caption <- paste0(
      "total = ", nrow(de),
      "; down = ", num_down,
      "; up = ", num_up
    )

    .plt <- make_volcano_plot(
      data = de,
      log_fc = avg_log2FC,
      log_pval = -log10(p_val_adj),
      direction = direction,
      gene_labels = tf_label
    ) + 
      labs(
        caption = caption,
        x = "log2(FC)",
        y = "-log10(adj. p-value)"
      ) + 
      scale_y_continuous(expand = expansion(mult = c(0.05, 0.1))) + 
      scale_colour_manual(
        values = c(up = "red", down = "blue", `n.s.` = "black"),
        guide = guide_legend(override.aes = list(alpha = 1, size = 3))
      ) + 
      theme_classic2() + 
      theme(legend.position = "bottom", legend.box.spacing = unit(1, "points"))

    ggsave(
      filename = paste0(subtype, "_TFs_volcano.png"),
      plot = .plt,
      width = 5,
      height = 5,
      units = "in",
      dpi = 600
    )
  }
)

# ------------------------------------------------------------------------------
# SOX4/SOX11 expression in MB subtypes (violin plot)

srat$mb$subtype <- factor(
  srat$mb$subtype,
  levels = c("SHH", "G3", "G4")
)

# prep for plotting
srat$mb <- PrepSCTFindMarkers(srat$mb)

# make violin plots
walk(
  .x = c("SOX4", "SOX11"),
  .f = \(gene) {
    .plt <- VlnPlot(
      srat$mb,
      features = gene,
      group.by = "subtype",
      assay = "SCT"
    ) + 
      theme_classic2()

    ggsave(
      filename = paste0(gene, "_violin.png"),
      plot = .plt,
      width = 6,
      height = 5,
      units = "in",
      dpi = 600
    )
  }
)

# ------------------------------------------------------------------------------
# heatmap of shared top regulons between UBCs and MB (code copied from
# 2024-09-29 notebook)

shared_regulons_dir <- "/.mounts/labs/pailab/private/icheong/CBL_scRNAseq/results/pyscenic/20240929/"
mat <- read.csv(
  file.path(shared_regulons_dir, "shared_regulons_mat.csv"),
  row.names = 1
) %>%
  as.matrix()

# colour function for heatmap
col_fun <- circlize::colorRamp2(
  breaks = c(0, max(mat)),
  hcl_palette = "Viridis"
)

# row/column title format
title_gp <- gpar(fontface = "bold")

# row/column names format
names_gp <- gpar(fontsize = 10)

# row labels
row_labels <- structure(
  rownames(mat) %>% str_replace_all(pattern = "_", replacement = " "),
  names = rownames(mat)
)

# main heatmap
hm <- Heatmap(
  matrix = mat,
  col = col_fun,
  name = "shared\nregulons",
  rect_gp = gpar(col = "white", lwd = 2),
  row_title = "UBC subclusters",
  row_title_side = "left",
  row_title_gp = title_gp,
  column_title = "MB subgroups",
  column_title_side = "bottom",
  column_title_gp = title_gp,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  row_labels = row_labels[rownames(mat)],
  row_names_side = "left",
  row_names_max_width = max_text_width(rownames(mat)),
  row_names_gp = names_gp,
  column_names_rot = 0,
  column_names_gp = names_gp,
  column_names_centered = TRUE,
  width = ncol(mat) * unit(0.5, "in"),
  height = nrow(mat) * unit(0.5, "in")
)
# need `draw()` to calculate heatmap width/height
pdf(NULL) # don't auto-generate "Rplots.pdf" file
hm <- draw(hm, merge_legend = TRUE, padding = unit(rep(0.5, 4), "in"))
dev.off()
# save heatmap to file
png(
  filename = "shared_regulons_heatmap.png",
  width = ComplexHeatmap:::width(hm),
  height = ComplexHeatmap:::height(hm),
  units = "mm",
  res = 600
)
# draw(hm, annotation_legend_list = lgd_list, merge_legend = TRUE)
draw(hm, merge_legend = TRUE)
dev.off()
