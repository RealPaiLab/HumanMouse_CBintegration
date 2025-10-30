# ==============================================================================
# These figures were generated for my committee meeting on 2024-07-24. This
# script should be run in the `scrnaseq_env` conda environment.
# ==============================================================================

library(tidyverse)
library(Seurat)
library(patchwork)
library(ggrepel)
library(Mfuzz)

# import functions
source("/.mounts/labs/pailab/private/icheong/CBL_scRNAseq/software/utilities/cell_labelling.R")
source("/.mounts/labs/pailab/private/icheong/CBL_scRNAseq/software/utilities/plotting.R")
source("/.mounts/labs/pailab/private/icheong/CBL_scRNAseq/software/utilities/propeller_helpers.R")

# >>> load Seurat objects >>>

# Seurat objects for import
srat_names <- c("rl_cca", "rl_rpca")
srat_qs <- c(
  # RL lineage CCA integration
  "/.mounts/labs/pailab/private/llau/results/integrated/20240524/25_pc_without_luo/25_pc_rl.qs",
  # RL lineage RPCA integration
  "/.mounts/labs/pailab/private/llau/results/integrated/20240606/rpca_20_analysis/rl/rl.qs"
)

# import Seurat objects
all_srat <- purrr::map(
  .x = srat_qs,
  .f = qs::qread
)
names(all_srat) <- srat_names
print(all_srat)

# from `cell_labelling.R`
all_srat$rl_cca <- label_rl_lineage_integration(all_srat$rl_cca)

# collapse UBC subclusters into one cluster
all_srat$rl_cca@meta.data <- mutate(
  all_srat$rl_cca[[]],
  broad_annot = case_when(
    str_starts(cell_type_annot, "UBC") ~ "UBC",
    .default = cell_type_annot
  )
)

# <<< load Seurat objects <<<

# >>> load gene expression trajectories >>>

# scaled expression of genes
gene_traj_df <- list()
for (lineage in c("ubc", "gc")) {
  lineage_dir <- file.path(
    "/.mounts/labs/pailab/private/icheong/CBL_scRNAseq/results/integrated/20240711",
    lineage
  )

  # scaled pseudobulk matrix
  pb_scaled <- readRDS(file.path(lineage_dir, "scaled_pseudobulk_eset.rds")) %>%
    exprs() %>%
    as.data.frame() %>%
    rownames_to_column("gene") %>%
    pivot_longer(
      cols = !gene,
      names_to = "bin",
      values_to = "scaled_expr"
    )

  # Mfuzz clustering membership scores
  mf_scores <- read_csv(file.path(lineage_dir, "mfuzz_membership_scores.csv")) %>%
    mutate(gene = paste(species, gene, sep = "_")) %>%
    dplyr::select(-species)

  # plot gene trajectories for each cluster
  gene_traj_df[[lineage]] <- left_join(
    x = pb_scaled,
    y = mf_scores,
    by = join_by(gene)
  )
}

# classification of genes as conserved/diverged
gene_class <- map(
  .x = c("ubc", "gc"),
  .f = \(lineage) {
    lineage_dir <- file.path(
      "/.mounts/labs/pailab/private/icheong/CBL_scRNAseq/results/integrated/20240711",
      lineage
    )
    df <- read_csv(file.path(lineage_dir, "mfuzz_clustering_pval.csv")) %>%
      mutate(lineage = lineage)
  }
) %>%
  do.call(rbind, args = .) %>%
  mutate(
    mf_class = case_when(
      pval < 0.05 & cluster_human != cluster_mouse ~ "diverged",
      pval > 0.5 & cluster_human == cluster_mouse ~ "conserved",
      .default = NA
    ),
    lineage = factor(lineage, levels = c("ubc", "gc"))
  )

# <<< load gene expression trajectories <<<

# colour palette for species
species_cols <- c("grey", scales::pal_hue()(1))

# set defaults for saving plots
hw_unit <- "in" # unit for height and width
dpi <- 600

# ------------------------------------------------------------------------------
# UMAPs of RL lineage CCA integration

rl_cluster <- "snn_res.0.4"
rl_cluster_cols <- pals::cols25(length(table(all_srat$rl_cca[[rl_cluster]])))
broad_annot_cols <- pals::trubetskoy(
  length(table(all_srat$rl_cca$broad_annot))
) %>%
  unname()
ubc_subcluster_cols <- pals::brewer.dark2(
  length(table(all_srat$rl_cca$ubc_subclusters))
)

# coloured by cluster
.plt <- DimPlot(
  all_srat$rl_cca,
  reduction = "umap",
  group.by = rl_cluster,
  cols = rl_cluster_cols,
  label = TRUE,
  repel = TRUE
) + 
  NoLegend() + 
  labs(title = NULL)
ggsave(
  file = "rl_cca_cluster.png",
  plot = .plt,
  width = 5.5,
  height = 5,
  units = hw_unit,
  dpi = dpi
)

# coloured by species
.plt <- DimPlot(
  all_srat$rl_cca,
  reduction = "umap",
  cells.highlight = WhichCells(all_srat$rl_cca, expression = species == "human"),
  sizes.highlight = 0.01,
  label = FALSE
) + 
  labs(title = NULL) + 
  scale_colour_manual(
    labels = c("mouse", "human"),
    values = species_cols
  )
ggsave(
  file = "rl_cca_species.png",
  plot = .plt,
  width = 6,
  height = 5,
  units = hw_unit,
  dpi = dpi
)

# coloured by broad cell type
.plt <- DimPlot(
  all_srat$rl_cca,
  reduction = "umap",
  group.by = "broad_annot",
  cols = broad_annot_cols,
  label = TRUE,
  repel = TRUE
) + 
  NoLegend() + 
  labs(title = NULL)
ggsave(
  file = "rl_cca_broad_annot.png",
  plot = .plt,
  width = 5.5,
  height = 5,
  units = hw_unit,
  dpi = dpi
)

# coloured by UBC subclusters
.plt <- DimPlot(
  all_srat$rl_cca,
  reduction = "umap",
  group.by = "ubc_subclusters",
  label = TRUE,
  repel = TRUE,
  cols = ubc_subcluster_cols,
  na.value = "grey"
) + 
  NoLegend() + 
  labs(title = NULL)
ggsave(
  file = "rl_ubc_subclusters.png",
  plot = .plt,
  width = 5.5,
  height = 5,
  units = hw_unit,
  dpi = dpi
)

# plot EOMES and LMX1A expression
walk(
  .x = c("EOMES", "LMX1A"),
  .f = \(gene) {
    .plt <- FeaturePlot(
      all_srat$rl_cca,
      features = gene,
      order = TRUE,
      min.cutoff = "q10",
      max.cutoff = "q90"
    )
    ggsave(
      filename = paste0(gene, ".png"),
      plot = .plt,
      width = 5,
      height = 4,
      units = hw_unit,
      dpi = dpi
    )
  }
)

# ------------------------------------------------------------------------------
# UMAPs of RL lineage RPCA integration

rl_cluster <- "snn_res.0.2"
rl_cluster_cols <- pals::cols25(length(table(all_srat$rl_rpca[[rl_cluster]])))

# coloured by cluster
.plt <- DimPlot(
  all_srat$rl_rpca,
  reduction = "umap",
  group.by = rl_cluster,
  cols = rl_cluster_cols,
  label = TRUE,
  repel = TRUE
) + 
  NoLegend() + 
  labs(title = NULL)
ggsave(
  file = "rl_rpca_cluster.png",
  plot = .plt,
  width = 5.5,
  height = 5,
  units = hw_unit,
  dpi = dpi
)

# coloured by species
.plt <- DimPlot(
  all_srat$rl_rpca,
  reduction = "umap",
  cells.highlight = WhichCells(all_srat$rl_rpca, expression = species == "human"),
  sizes.highlight = 0.01,
  label = FALSE
) + 
  labs(title = NULL) + 
  scale_colour_manual(
    labels = c("mouse", "human"),
    values = species_cols
  )
ggsave(
  file = "rl_rpca_species.png",
  plot = .plt,
  width = 6,
  height = 5,
  units = hw_unit,
  dpi = dpi
)

# ------------------------------------------------------------------------------
# differential abundance (cell type proportions)

scpt <- list()

# function for plotting scProportionTest results
plot_scpt <- function(
  scpt_df,
  order_clusters = TRUE # order by log2FD value
) {
  if (order_clusters) {
    scpt_df$clusters <- fct_reorder(factor(scpt_df$clusters), scpt_df$obs_log2FD)
  }

  # plot log2FD for each cluster
  p <- ggplot(
    data = scpt_df,
    mapping = aes(x = clusters, y = obs_log2FD)
  ) + 
    geom_pointrange(aes(
      ymin = boot_CI_2.5,
      ymax = boot_CI_97.5,
      colour = significance
    )) + 
    geom_hline(yintercept = 0, linetype = "solid") + 
    geom_hline(yintercept = c(-1, 1), linetype = "dashed", colour = "grey50") + 
    scale_colour_manual(values = c("#007894", "grey")) + 
    coord_flip() + 
    theme_light() + 
    theme(
      axis.text = element_text(colour = "black")
    )
  
  return(p)
}

# cell type annotations for each cluster
cluster_annot <- read_csv("/.mounts/labs/pailab/private/llau/results/integrated/20240618/clusters_cell_types.csv")

broad_annot <- cluster_annot %>%
  mutate(
    final_clusters = case_when(
      str_detect(Cluster, pattern = "Cluster_") ~ str_remove(Cluster, "^Cluster_"),
      str_detect(Cluster, pattern = "UBC_") ~ "2",
      .default = Cluster
    ),
    .before = 1
  ) %>%
  dplyr::rename(broad_annot = broad_cell_type) %>%
  dplyr::select(final_clusters, broad_annot) %>%
  dplyr::distinct()

ubcsub_annot <- cluster_annot %>%
  mutate(final_clusters = str_remove(Cluster, "^Cluster_")) %>%
  dplyr::rename(cell_type_annot = broad_w_ubc_subtypes) %>%
  dplyr::select(final_clusters, cell_type_annot)

# >>> RL clusters (no UBC subclusters) >>>

# load scProportionTest results and add cell type annotations
scpt$broad <- read_csv("/.mounts/labs/pailab/private/llau/results/integrated/20240718/proportion_test/cca/cca_rl_clusters_all_timepoints_proportion_test.csv") %>%
  mutate(
    clusters = as.character(clusters),
    significance = case_when(
      FDR < 0.05 & abs(obs_log2FD) > 1 ~ "FDR < 0.05 & |log2FD| > 1",
      .default = "not sig."
    )
  ) %>%
  left_join(
    x = .,
    y = broad_annot,
    by = join_by(clusters == final_clusters)
  ) %>% 
  mutate(
    clusters = paste0(clusters, " (", broad_annot, ")")
  )

# plot scProportionTest
.plt <- plot_scpt(scpt_df = scpt$broad) + 
  theme(legend.position = "bottom")
ggsave(
  filename = "scpt_broad_clusters.png",
  plot = .plt,
  width = 5,
  height = 6,
  units = hw_unit,
  dpi = dpi
)

# load propeller results and filter for UBCs
pplr_broad_prop <- read_csv(
  "/.mounts/labs/pailab/private/icheong/CBL_scRNAseq/results/integrated/20240617/cca_proportions_without_ubc_subclusters.csv"
) %>%
  mutate(species = factor(species, levels = c("mouse", "human")))
pplr_broad_res <- read_csv(
  "/.mounts/labs/pailab/private/icheong/CBL_scRNAseq/results/integrated/20240617/cca_results_without_ubc_subclusters.csv"
)

# plot propeller results
broad_fdrs <- get_fdrs_for_plotting(pplr_broad_res)
.plt <- boxplot_props_per_cluster(
  props = pplr_broad_prop %>% filter(clusters == "2"),
  facet_var = clusters,
  facet_ncol = 3,
  fdrs = broad_fdrs %>% filter(clusters == "2")
) + 
  theme(axis.text = element_text(colour = "black"))
ggsave(
  "propeller_broad_clusters.png",
  plot = .plt,
  width = 4,
  height = 5,
  units = hw_unit,
  dpi = dpi
)

# <<< RL clusters (no UBC subclusters) <<<

# >>> RL clusters (with UBC subclusters) >>>

# bar plot with proportion of cells in each UBC subcluster
rl_cca <- all_srat$rl_cca
rl_cca@meta.data <- mutate(
  rl_cca[[]],
  ubc_subclusters = case_when(
    is.na(ubc_subclusters) ~ "non-UBC",
    .default = paste0("UBC_", ubc_subclusters)
  ),
  ubc_subclusters = fct_relevel(factor(ubc_subclusters), "non-UBC", after = Inf)
)
.plt <- cluster_barplot(
  object = rl_cca,
  split.by = "ubc_subclusters",
  group.by = "species",
  position = "fill"
) + 
  scale_fill_manual(values = c(ubc_subcluster_cols, "grey50")) + 
  coord_cartesian(ylim = c(0.7, 1)) + 
  theme(legend.title = element_blank())
ggsave(
  filename = "ubc_subclusters_bar.png",
  plot = .plt,
  width = 3.5,
  height = 5,
  units = hw_unit,
  dpi = dpi
)

# load scProportionTest results and add cell type annotations
scpt$ubcsub <- read_csv("/.mounts/labs/pailab/private/llau/results/integrated/20240613/cca/cca_rl_proportion_test.csv") %>%
  mutate(
    clusters = str_remove(clusters, "^Cluster_"),
    significance = case_when(
      FDR < 0.05 & abs(obs_log2FD) > 1 ~ "FDR < 0.05 & |log2FD| > 1",
      .default = "not sig."
    )
  ) %>%
  left_join(
    x = .,
    y = ubcsub_annot,
    by = join_by(clusters == final_clusters)
  ) %>% 
  mutate(
    clusters = case_when(
      !str_detect(clusters, "^UBC_") ~ paste0(clusters, " (", cell_type_annot, ")"),
      .default = clusters
    )
  )

# plot scProportionTest
.plt <- plot_scpt(scpt_df = scpt$ubcsub) + 
  theme(legend.position = "bottom")
ggsave(
  filename = "scpt_ubc_subclusters.png",
  plot = .plt,
  width = 5,
  height = 6,
  units = hw_unit,
  dpi = dpi
)

# load propeller results and filter for UBCs
pplr_ubcsub_prop <- read_csv(
  "/.mounts/labs/pailab/private/icheong/CBL_scRNAseq/results/integrated/20240617/cca_proportions_with_ubc_subclusters.csv"
) %>%
  mutate(species = factor(species, levels = c("mouse", "human")))
pplr_ubcsub_res <- read_csv(
  "/.mounts/labs/pailab/private/icheong/CBL_scRNAseq/results/integrated/20240617/cca_results_with_ubc_subclusters.csv"
)

# plot propeller results
ubcsub_fdrs <- get_fdrs_for_plotting(pplr_ubcsub_res)
.plt <- boxplot_props_per_cluster(
  props = pplr_ubcsub_prop %>% filter(str_detect(clusters, "^UBC_")),
  facet_var = clusters,
  facet_ncol = 3,
  fdrs = ubcsub_fdrs %>% filter(str_detect(clusters, "UBC_"))
) + 
  theme(axis.text = element_text(colour = "black"))
ggsave(
  "propeller_ubc_subclusters.png",
  plot = .plt,
  width = 7,
  height = 8,
  units = hw_unit,
  dpi = dpi
)

# <<< RL clusters (with UBC subclusters) <<<

# ------------------------------------------------------------------------------
# differential gene expression

diff_exp <- read_csv("/.mounts/labs/pailab/private/llau/results/integrated/20240715/all_tested_genes.csv") %>%
  # filter for UBC_1
  filter(cluster == "UBC_1")

# clean up column names for easier access
colnames(diff_exp) <- str_remove(colnames(diff_exp), "full_cerebellum_")

# add gene labels for plotting
top_aldinger_features <- diff_exp %>%
  arrange(desc(abs(Aldinger_human_avg_log2FC))) %>%
  head(n = 25) %>%
  pull(feature)
top_sepp_features <- diff_exp %>%
  arrange(desc(abs(Sepp_human_avg_log2FC))) %>%
  head(n = 25) %>%
  pull(feature)
diff_exp <- mutate(
  diff_exp,
  labels = case_when(
    feature %in% c(top_aldinger_features, top_sepp_features) ~ feature,
    .default = ""
  ),
  Aldinger_direction = case_when(
    Aldinger_human_avg_log2FC > 0 ~ "up",
    Aldinger_human_avg_log2FC < 0 ~ "down"
  ),
  Sepp_direction = case_when(
    Sepp_human_avg_log2FC > 0 ~ "up",
    Sepp_human_avg_log2FC < 0 ~ "down"
  )
)

# make volcano plot for Aldinger dataset
.plt <- make_volcano(
  data = diff_exp,
  log_fc = Aldinger_human_avg_log2FC,
  log_pval = -log10(Aldinger_human_p_val_adj),
  direction = Aldinger_direction,
  gene_labels = labels,
  log_pval_thresh = -log10(0.05),
  label_num_genes = TRUE
) + 
  labs(x = "log2(FC)", y = "-log10(adj. p-value)") + 
  theme_light() + 
  theme(
    axis.text = element_text(colour = "black"),
    legend.position = "bottom"
  )
ggsave(
  filename = "volcano_aldinger.png",
  plot = .plt,
  width = 5,
  height = 6,
  units = hw_unit,
  dpi = dpi
)

# make volcano plot for Sepp dataset
.plt <- make_volcano(
  data = diff_exp,
  log_fc = Sepp_human_avg_log2FC,
  log_pval = -log10(Sepp_human_p_val_adj),
  direction = Sepp_direction,
  gene_labels = labels,
  log_pval_thresh = -log10(0.05),
  label_num_genes = TRUE
) + 
  labs(x = "log2(FC)", y = "-log10(adj. p-value)") + 
  theme_light() + 
  theme(
    axis.text = element_text(colour = "black"),
    legend.position = "bottom"
  )
ggsave(
  filename = "volcano_sepp.png",
  plot = .plt,
  width = 5,
  height = 6,
  units = hw_unit,
  dpi = dpi
)

# ------------------------------------------------------------------------------
# gene expression trajectories

# function to plot gene trajectories
plot_gene_traj <- function(
  df,
  cluster_score = NULL
) {
  p <- ggplot(
    data = df,
    mapping = aes(x = bin, y = scaled_expr, group = gene)
  )

  if (!is.null(cluster_score)) {
    p <- p + geom_line(aes(colour = .data[[cluster_score]]))
  } else {
    p <- p + geom_line()
  }

  p <- p + 
    labs(x = "pseudotime bin", y = "scaled expression") + 
    theme_classic() + 
    theme(
      axis.text = element_text(colour = "black"),
      axis.ticks = element_line(colour = "black")
    )
  return(p)
}

# plot gene trajectories for the UBC and GC lineages
walk2(
  .x = gene_traj_df,
  .y = names(gene_traj_df),
  .f = \(gene_traj_df, lineage) {
    .plt_list <- map(
      .x = c(1:4),
      .f = \(cluster) {
        plot_title <- paste("cluster", cluster)
        cluster <- paste0("mf_mem_", cluster)
        gene_traj_df <- gene_traj_df %>%
          # keep genes with a membership score > 0.5 in that column
          dplyr::filter(.data[[cluster]] > 0.5)

        # set factor level so genes with higher scores get plotted on top
        gene_lvl <- gene_traj_df[, c("gene", cluster)] %>%
          distinct() %>%
          arrange(.data[[cluster]]) %>%
          pull(gene)
        gene_traj_df <- mutate(
          gene_traj_df,
          gene = factor(gene, levels = gene_lvl)
        )

        # plot gene trajectories
        p <- plot_gene_traj(gene_traj_df, cluster_score = cluster) + 
          labs(title = plot_title, colour = "membership\nscore")

        return(p)
      }
    )

    # combine plots
    .plt <- wrap_plots(.plt_list, ncol = 2, guides = "collect") & 
      scale_colour_viridis_c() & 
      theme(plot.title = element_text(face = "bold", hjust = 0.5))

    ggsave(
      filename = paste0(lineage, "_gene_traj_clusters.png"),
      plot = .plt,
      width = 7,
      height = 6,
      units = hw_unit,
      dpi = dpi
    )
  }
)

# plot number of diverged/conserved genes in UBC/GC lineages
gene_class_summary <- gene_class %>%
  group_by(lineage, mf_class) %>%
  summarise(n = n()) %>%
  ungroup(mf_class) %>%
  mutate(prop = n / sum(n))

.plt <- ggplot(
  data = gene_class_summary,
  mapping = aes(x = lineage, y = n, fill = mf_class, label = n)
) + 
  geom_col(width = 0.6) + 
  geom_text(position = position_stack(vjust = 0.5)) + 
  labs(x = "lineage", y = "number of genes") + 
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) + 
  scale_fill_hue(l = 80, na.value = "grey75") + 
  theme_light() + 
  theme(
    axis.text = element_text(colour = "black"),
    legend.position = "bottom"
  )
ggsave(
  filename = "gene_traj_class_bar.png",
  plot = .plt,
  width = 4,
  height = 5,
  units = hw_unit,
  dpi = dpi
)

# plot PKNOX2 expression
.plt <- ggplot(
  data = gene_traj_df$gc %>% 
    filter(str_detect(gene, pattern = "PKNOX2")) %>%
    separate_wider_delim(gene, delim = "_", names = c("species", "gene")),
  mapping = aes(x = bin, y = scaled_expr, group = species, colour = species)
) + 
  geom_line() + 
  labs(x = "pseudotime bin", y = "scaled expression") + 
  scale_colour_manual(values = rev(species_cols)) + 
  theme_classic() + 
  theme(
    axis.text = element_text(colour = "black"),
    axis.ticks = element_line(colour = "black")
  )
ggsave(
  filename = "PKNOX2_traj.png",
  plot = .plt,
  width = 4,
  height = 3,
  units = hw_unit,
  dpi = dpi
)

.plt <- FeaturePlot(
  all_srat$rl_cca,
  features = "PKNOX2",
  order = TRUE,
  min.cutoff = "q10",
  max.cutoff = "q90",
  split.by = "species"
)
ggsave(
  filename = "PKNOX2.png",
  plot = .plt,
  width = 8,
  height = 4,
  units = hw_unit,
  dpi = dpi
)
