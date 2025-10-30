# ==============================================================================
# Figure 2 of the manuscript. This script should be run using the scrnaseq_env
# conda environment.
# ==============================================================================

library(tidyverse)
library(patchwork)
library(Seurat)
library(ggalluvial)

source("./utils.R")
source("../../software/utilities/cell_labelling.R")
source("../../software/utilities/plotting.R")
source("../../software/utilities/propeller_helpers.R")

out_dir <- sprintf("/.mounts/labs/pailab/private/projects/HumanMouseUBC/figures/Fig2_human_mouse_integ_%s", format(Sys.Date(), "%Y%m%d"))
my_pals <- get_custom_pals()

plt <- list()

# ------------------------------------------------------------------------------
# UMAPs of integrated human/mouse RL lineage cells

# load Seurat object for integrated dataset
srat_qs <- get_srat_paths()
rl_srat <- load_srat(srat_qs["rl"]) %>%
  pluck("rl")

# get cell type annotations (from `cell_labelling.R`)
rl_srat <- label_rl_lineage_integration(rl_srat)

# collapse UBC subclusters into one cluster
rl_srat@meta.data <- mutate(
  rl_srat[[]],
  broad_annot = case_when(
    str_starts(cell_type_annot, "UBC") ~ "UBC",
    .default = cell_type_annot
  )
)

# set column for clusters
rl_cluster <- "snn_res.0.4"

# UMAP coloured by species
plt$umap_rl_species <- DimPlot(
  rl_srat,
  reduction = "umap",
  cells.highlight = WhichCells(rl_srat, expression = species == "human"),
  sizes.highlight = 0.01,
  label = FALSE
) + 
  labs(title = "Species", x = "UMAP 1", y = "UMAP 2") + 
  scale_colour_manual(
    labels = c("mouse", "human"),
    values = my_pals$species,
    guide = guide_legend(override.aes = list(size = 3), nrow = 1)
  ) + 
  theme_classic2() + 
  theme(
    legend.position = "bottom",
    legend.justification.bottom = "top"
  )
ggsave(
  filename = "integ_rl_species.pdf",
  plot = plt$umap_rl_species,
  path = out_dir,
  width = 3.5,
  height = 3.5,
  units = "in"
)

# UMAPs coloured by cluster and cell type
plt[c("umap_rl_cell_type", "umap_rl_cluster")] <- map(
  .x = c("broad_annot", rl_cluster),
  .f = \(group_by) {
    if (str_detect(group_by, "snn_res")) {
      cols <- my_pals$rl_integ_clust
      filename <- "integ_rl_cluster.pdf"
      title <- "Cluster"
    } else if (group_by == "broad_annot") {
      cols <- my_pals$rl_integ_annot
      filename <- "integ_rl_broad_annot.pdf"
      title <- "Rhombic lip"
    }

    .plt <- DimPlot(
      rl_srat,
      reduction = "umap",
      group.by = group_by,
      cols = cols,
      label = TRUE,
      label.size = 3,
      repel = TRUE
    ) + 
      labs(title = title, x = "UMAP 1", y = "UMAP 2") + 
      guides(colour = guide_legend(
        override.aes = list(size = 3),
        nrow = 3
      )) + 
      theme_classic2() + 
      theme(legend.position = "bottom")

    # hacky way of setting the seed for the UMAP labels
    .plt[[1]]$layers[[2]]$geom_params$seed <- 42

    ggsave(
      filename = filename,
      plot = .plt,
      path = out_dir,
      width = 3.5,
      height = 4,
      units = "in"
    )

    return(.plt)
  }
)

# expression of EOMES on UMAP
plt$umap_rl_eomes <- FeaturePlot(
  rl_srat,
  features = "EOMES",
  order = TRUE,
  min.cutoff = "q10",
  max.cutoff = "q90",
  combine = TRUE
) + 
  labs(x = "UMAP 1", y = "UMAP 2") + 
  scale_colour_gradient(
    low = "lightgrey",
    high = "blue",
    breaks = \(x) {x}, # use min/max expression as breaks
    labels = c("low", "high")
  ) + 
  theme_classic2() + 
  theme(
    legend.ticks = element_blank(),
    legend.position = "bottom",
    legend.justification.bottom = "top"
  )
ggsave(
  filename = "integ_rl_EOMES.pdf",
  plot = plt$umap_rl_eomes,
  path = out_dir,
  width = 3.5,
  height = 3.75,
  units = "in"
)

# ------------------------------------------------------------------------------
# UMAPs of integrated human/mouse UBCs

ubc_srat <- load_srat(srat_qs["ubc"]) %>%
  pluck("ubc")

ubc_srat$subclusters <- paste0("i", str_remove(
  string = ubc_srat$subclusters,
  pattern = "_"
))

# UMAP coloured by subclusters
plt$umap_ubc_cluster <- DimPlot(
  ubc_srat,
  reduction = "umap",
  group.by = "subclusters",
  cols = my_pals$ubc_integ_clust2,
  label = TRUE,
  label.size = 3,
  repel = TRUE
) + 
  labs(title = "UBC clusters", x = "UMAP 1", y = "UMAP 2") + 
  theme_classic2() + 
  theme(legend.position = "none")

# hacky way of setting the seed for the UMAP labels
plt$umap_ubc_cluster[[1]]$layers[[2]]$geom_params$seed <- 42

ggsave(
  filename = "integ_ubc_cluster.pdf",
  plot = plt$umap_ubc_cluster,
  path = out_dir,
  width = 3.5,
  height = 3.5,
  units = "in"
)

# UMAP coloured by species
plt$umap_ubc_species <- DimPlot(
  ubc_srat,
  reduction = "umap",
  pt.size = 0.2,
  cells.highlight = WhichCells(ubc_srat, expression = species == "human"),
  sizes.highlight = 0.2,
  label = FALSE
) + 
  labs(title = "Species", x = "UMAP 1", y = "UMAP 2") + 
  scale_colour_manual(
    labels = c("mouse", "human"),
    values = my_pals$species,
    guide = guide_legend(override.aes = list(size = 3), nrow = 1)
  ) + 
  theme_classic2() + 
  theme(
    legend.position = "inside",
    legend.position.inside = c(0.01, 1),
    legend.justification.inside = c(0, 1)
  )
ggsave(
  filename = "integ_ubc_species.pdf",
  plot = plt$umap_ubc_species,
  path = out_dir,
  width = 3.5,
  height = 3.5,
  units = "in"
)

# ------------------------------------------------------------------------------
# UBC proportions

# set species factor levels
ubc_srat$species <- factor(ubc_srat$species, levels = c("mouse", "human"))

# proportion of UBC subclusters in human vs mouse
plt$prop_bar <- cluster_barplot(
  object = ubc_srat,
  split.by = "subclusters",
  group.by = "species",
  position = "fill"
) + 
  labs(fill = "subcluster") + 
  scale_fill_manual(values = my_pals$ubc_integ_clust2) + 
  theme_classic2()
ggsave(
  filename = "ubc_cluster_prop_bar.pdf",
  plot = plt$prop_bar,
  path = out_dir,
  width = 3,
  height = 3.5,
  units = "in"
)

# >>> scProportionTest >>>

# load scProportionTest results
scpt_res <- read_csv("/.mounts/labs/pailab/private/llau/results/integrated/20240718/proportion_test/cca/cca_ubc_clusters_proportion_test.csv")

# add significance column and set plot order
log2FD_threshold <- log2(2)
FDR_threshold <- 0.05
scpt_res <- mutate(
  scpt_res,
  significance = case_when(
    FDR < FDR_threshold & obs_log2FD > log2FD_threshold ~ "human",
    FDR < FDR_threshold & obs_log2FD < -log2FD_threshold ~ "mouse",
    .default = "n.s."
  ),
  significance = factor(significance, levels = c("mouse", "n.s.", "human")),
  clusters = paste0("i", str_remove(clusters, "_")),
  clusters = fct_reorder(clusters, obs_log2FD, .desc = FALSE)
)

# plot scProportionTest results
plt$scpt <- ggplot(
  data = scpt_res,
  mapping = aes(
    x = obs_log2FD,
    y = clusters,
    xmin = boot_CI_2.5,
    xmax = boot_CI_97.5,
    colour = significance
  )
) + 
  geom_pointrange(size = 0.25) + 
  geom_vline(
    xintercept = c(log2FD_threshold, -log2FD_threshold),
    lty = "dashed"
  ) + 
  geom_vline(xintercept = 0, lty = "solid") + 
  labs(
    title = "scProportionTest for UBC clusters",
    x = "observed log2(FD)",
    y = "UBC cluster",
    colour = "enriched"
  ) + 
  scale_colour_manual(values = c(my_pals$species[1], "black", my_pals$species[2])) + 
  theme_classic2() + 
  theme(legend.position = "bottom", plot.title = element_text(vjust = 2))
ggsave(
  filename = "ubc_permutation_plot.pdf",
  plot = plt$scpt,
  path = out_dir,
  width = 4,
  height = 3,
  units = "in"
)

# <<< scProportionTest <<<

# >>> propeller >>>

# set same cluster plotting order as scProportionTest
cluster_levels <- sort(levels(scpt_res$clusters))

# load `propeller` sample proportions
propeller_props <- read_csv("/.mounts/labs/pailab/private/icheong/CBL_scRNAseq/results/integrated/20240606/four_dataset_cca_ubc_proportions.csv") %>%
  mutate(
    dataset_name = str_replace(
      string = dataset_name,
      pattern = "_([:alpha:]+)",
      replacement = " (\\1)"
    ),
    species = factor(species, levels = c("mouse", "human")),
    clusters = paste0("UBC_", clusters) %>% factor(levels = cluster_levels)
  )

# load propeller
propeller_res <- read_csv("/.mounts/labs/pailab/private/icheong/CBL_scRNAseq/results/integrated/20240606/four_dataset_cca_ubc_results.csv")
fdrs <- get_fdrs_for_plotting(propeller_res) %>%
  mutate(
    clusters = paste0("UBC_", clusters) %>% factor(levels = cluster_levels)
  )

# plot `propeller` results for UBCs
plt$propeller <- boxplot_props_per_cluster(
  props = propeller_props,
  facet_var = clusters,
  facet_ncol = 6,
  fdrs = fdrs,
  fdr_label_size = 8.8
) + 
  labs(y = "proportion of all\nRL lineage cells", shape = "dataset") + 
  theme_classic2()
ggsave(
  filename = "ubc_propeller_plot.pdf",
  plot = plt$propeller,
  path = out_dir,
  width = 14,
  height = 4,
  units = "in"
)

# <<< propeller <<<

# ------------------------------------------------------------------------------
# map clusters back to human-only integration (Aldinger/Sepp UBCs)

quang_ubc_srat <- load_srat(srat_qs["quang_ubc"]) %>%
  pluck("quang_ubc")

# add Quang's clusters to metadata
ubc_srat$human_only_clusters <- quang_ubc_srat[[]][rownames(ubc_srat[[]]), "SCT_snn_res.0.5"]

# make bar plots
plt$mapping <- cluster_barplot(
  object = subset(ubc_srat, subset = species == "human"), # plot human cells only
  split.by = "human_only_clusters",
  group.by = "subclusters",
  position = "stack",
  width = 0.75
) + 
  labs(
    title = "Human UBCs in the integrated clusters\ncompared to human-only clustering",
    fill = "original human-only\nclusters"
  ) + 
  # scale_fill_manual(values = my_pals$ubc_integ_clust) + 
  theme_classic2()
ggsave(
  plot = plt$mapping,
  filename = "ubc_cluster_mappings_bar.pdf",
  path = out_dir,
  width = 6,
  height = 4,
  units = "in"
)

# ------------------------------------------------------------------------------
# volcano plots

# load differential gene expression results (all UBC clusters, both datasets)
ubc_diff_expr <- read_csv("/.mounts/labs/pailab/private/llau/results/integrated/20240715/all_tested_genes.csv") %>%
  rename_with(.fn = str_remove, pattern = "full_cerebellum_")

# load top markers in each cluster across both datasets (Aldinger/Sepp)
ranked_markers <- readRDS("/.mounts/labs/pailab/private/icheong/CBL_scRNAseq/results/integrated/20241029/cluster_marker_ranking.rds")

# volcano plots for UBC_1 and UBC_2
plt$volcano <- map(
  .x = paste0("UBC_", c(1:2)),
  .f = \(clust) {
    # top 10 features for labelling
    features <- ranked_markers %>%
      pluck(clust) %>%
      slice_min(combined_rank, n = 10, with_ties = TRUE) %>%
      pull(feature)

    # volcano plot colours
    volcano_cols <- c(
      "down" = my_pals$rl_integ_annot[["UBC"]],
      "n.s." = "black",
      "up" = my_pals$ubc_integ_clust[[clust]]
    )

    p <- list()
    for (dataset in c("Aldinger", "Sepp")) {
      # which column to use (based on dataset)
      logfc_column <- paste0(str_to_title(dataset), "_human_avg_log2FC")
      pval_column <- paste0(str_to_title(dataset), "_human_p_val")

      # filter for cluster and dataset to be plotted
      diff_expr <- ubc_diff_expr %>%
        filter(cluster == clust) %>%
        mutate(
          direction = case_when(
            rlang::parse_expr(
              paste0(logfc_column, " < 0 & ", pval_column, " < 0.05 ~ 'down'")
            ),
            rlang::parse_expr(
              paste0(logfc_column, " > 0 & ", pval_column, " < 0.05 ~ 'up'")
            ),
            .default = "n.s."
          ),
          gene_label = case_when(
            feature %in% features ~ feature,
            .default = ""
          ),
        )

      # label number of genes
      num_up <- nrow(filter(diff_expr, direction == "up"))
      num_down <- nrow(filter(diff_expr, direction == "down"))
      caption <- paste0(
        "total = ", nrow(diff_expr),
        "; down = ", num_down,
        "; up = ", num_up
      )

      # volcano plot
      p[[dataset]] <- make_volcano(
        data = diff_expr,
        log_fc = .data[[logfc_column]],
        log_pval = -log10(.data[[pval_column]]),
        direction = direction,
        gene_label = NULL,
        label_num_genes = FALSE
      ) + 
        ggrepel::geom_text_repel(
          aes(label = gene_label, colour = direction),
          size = 3.5,
          max.overlaps = Inf,
          show.legend = FALSE,
          seed = 42,
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
        ) + 
        labs(
          title = str_to_title(dataset),
          caption = caption,
          x = "log2(FC)",
          y = "-log10(p-value)"
        ) + 
        scale_colour_manual(values = volcano_cols) +  
        guides(colour = guide_legend(
          title = clust,
          override.aes = list(alpha = 1)
        )) + 
        theme_classic2()
    }

    # combine plots
    p <- wrap_plots(free(p$Aldinger), free(p$Sepp), ncol = 2, guides = "collect")
    ggsave(
      filename = paste0("volcano_", clust, ".pdf"),
      plot = p,
      path = out_dir,
      width = 8,
      height = 4,
      units = "in",
      dpi = 600
    )

    return(p)
  }
) %>%
  set_names(paste0("UBC_", c(1:2)))

# ------------------------------------------------------------------------------
# pathway enrichment analysis

# load pathway enrichment results for all UBC subclusters
ubc_pathways <- map(
  .x = paste0("UBC_", c(1:2)),
  .f = \(clust) {
    ubc_pathways <- read_csv(file.path(
      "/.mounts/labs/pailab/private/icheong/CBL_scRNAseq/results/integrated/20241029",
      clust,
      "enr_pathways.csv"
    )) %>%
      # capitalize term names
      mutate(term_name = str_to_sentence(term_name))

    return(ubc_pathways)
  }
) %>%
  setNames(nm = paste0("UBC_", c(1:2)))

# enriched pathways for all UBC subclusters
plt$pathway_bar <- map(
  .x = names(ubc_pathways),
  .f = \(clust) {
    # subset pathways for the subcluster
    pway <- ubc_pathways[[clust]]

    # number of significant pathways
    n_sig <- nrow(filter(pway, p_value < 0.05))

    # total number of enriched pathways
    caption <- paste0("# enriched pathways/total = ", n_sig, "/", nrow(pway))

    # how many pathways to show?
    if (n_sig == 0) {
      # if no significant pathways, just show top 10
      n_plot <- 10
    } else {
      # show up to 10 pathways
      n_plot <- min(n_sig, 10)
    }

    # get pathways to plot
    pway <- slice_min(
      pway,
      order_by = p_value,
      n = n_plot
    )

    # set font size and width of terms
    term_font_size <- rel(0.8)
    lab_width <- 30

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
        title = sprintf("Enriched pathways in %s", clust),
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
      filename = sprintf("enr_pathways_%s.pdf", clust),
      plot = .plt,
      path = out_dir,
      width = 4,
      height = 3.5,
      units = "in",
      dpi = 600
    )

    return(.plt)
  }
) %>%
  setNames(nm = names(ubc_pathways))

# ------------------------------------------------------------------------------
# combine plots

layout <- c(
  area(t = 1, b = 2, l = 1, r = 4),
  area(t = 1, b = 2, l = 5, r = 8),
  area(t = 1, b = 2, l = 9, r = 12),
  area(t = 3, b = 4, l = 1, r = 4),
  area(t = 3, b = 4, l = 5, r = 7),
  area(t = 3, b = 4, l = 8, r = 12),
  area(t = 5, b = 6, l = 1, r = 8),
  area(t = 7, b = 8, l = 1, r = 7),
  area(t = 7, b = 8, l = 8, r = 11),
  area(t = 9, b = 10, l = 1, r = 7),
  area(t = 9, b = 10, l = 8, r = 11)
)

.plt <- wrap_plots(
  plt$umap_rl_cell_type,
  plt$umap_rl_species,
  plt$umap_rl_eomes,
  plt$umap_ubc_cluster,
  free(plt$prop_bar),
  free(plt$scpt),
  free(plt$mapping),
  free(plt$volcano$UBC_1 & labs(caption = NULL)),
  free(plt$pathway_bar$UBC_1 + labs(caption = NULL)),
  free(plt$volcano$UBC_2 & labs(caption = NULL)),
  free(plt$pathway_bar$UBC_2 + labs(caption = NULL)),
  design = layout
)
walk(
  .x = c("png", "pdf"),
  .f = \(ext) {
    ggsave(
      filename = sprintf("combined_plots.%s", ext),
      plot = .plt,
      path = out_dir,
      width = 12,
      height = 18,
      units = "in",
      dpi = 600
    )
  }
)
