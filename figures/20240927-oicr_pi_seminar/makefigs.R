# ==============================================================================
# These figures were generated for the OICR PI seminar on 2024-09-24. This
# script should be run in the `scrnaseq_env` conda environment.
# ==============================================================================

library(tidyverse)
library(Seurat)
library(patchwork)

# import functions
source("/.mounts/labs/pailab/private/icheong/CBL_scRNAseq/software/utilities/cell_labelling.R")
source("/.mounts/labs/pailab/private/icheong/CBL_scRNAseq/software/utilities/plotting.R")
source("/.mounts/labs/pailab/private/icheong/CBL_scRNAseq/thesis/utils.R")

# colour palettes (function from `utils.R`)
my_pals <- get_custom_pals()

# add medulloblastoma subtypes to colour palettes
my_pals$mb_subtype <- pals::brewer.set2(5)[3:5]

# >>> load Seurat objects >>>

# `get_srat_paths` and `load_srat` from `utils.R`
srat_qs <- get_srat_paths()
srat <- load_srat(srat_qs[c("rl", "mb")])

# get cell type annotations (from `cell_labelling.R`)
srat$rl <- label_rl_lineage_integration(srat$rl)

# collapse UBC subclusters into one cluster
srat$rl@meta.data <- mutate(
  srat$rl[[]],
  broad_annot = case_when(
    str_starts(cell_type_annot, "UBC") ~ "UBC",
    .default = cell_type_annot
  ),
  broad_annot = factor(
    broad_annot,
    levels = c("RL", "UBC", "GCP", "GC", "oligodendrocyte/OPC", "microglia", "endothelial")
  ),
  dataset_name = str_remove(
    string = dataset_name,
    pattern = "full_cerebellum_"
  ),
  species = factor(species, levels = c("mouse", "human"))
)

# normalize the RNA assay (for dot plot)
DefaultAssay(srat$rl) <- "RNA"
srat$rl <- NormalizeData(srat$rl)

# set levels for medulloblastoma subtypes
mb_subtypes <- c("SHH", "G3", "G4")
srat$mb$subtype <- factor(srat$mb$subtype, levels = mb_subtypes)

# <<< load Seurat objects <<<

# load list of transcription factors
tfs <- read.csv(
  "/.mounts/labs/pailab/src/transcription-factors/human-tfs-lambert/full_database_v1.01.csv",
  row.names = 1,
  check.names = FALSE
) %>%
  dplyr::filter(`Is TF?` == "Yes") %>%
  pull(`HGNC symbol`)

# ------------------------------------------------------------------------------
# RL lineage UMAPs; code taken from `integ_rl.qmd` in thesis figures directory

# UMAP coloured by species
.plt <- DimPlot(
  srat$rl,
  reduction = "umap",
  cells.highlight = WhichCells(srat$rl, expression = species == "human"),
  sizes.highlight = 0.01,
  label = FALSE
) + 
  labs(title = NULL) + 
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
  filename = "rl_umap_species.png",
  plot = .plt,
  width = 4,
  height = 4,
  units = "in",
  dpi = 600
)

# UMAPs coloured by cluster and cell type
walk(
  .x = c("broad_annot", "snn_res.0.4"),
  .f = \(group_by) {
    if (str_detect(group_by, "snn_res")) {
      cols <- my_pals$rl_integ_clust
      legend_pos <- "none"
      filename <- "rl_umap_cluster.png"
      width <- 4
    } else if (group_by == "broad_annot") {
      cols <- my_pals$rl_integ_annot
      legend_pos <- "right"
      filename <- "rl_umap_cell_type.png"
      width <- 5.25
    }

    .plt <- DimPlot(
      srat$rl,
      reduction = "umap",
      group.by = group_by,
      cols = cols,
      label = TRUE,
      label.size = 3,
      repel = TRUE
    ) + 
      labs(title = NULL) + 
      guides(colour = guide_legend(
        override.aes = list(size = 3)
        # nrow = 3
      )) + 
      theme_classic2() + 
      theme(legend.position = legend_pos)

    # hacky way of setting the seed for the UMAP labels
    .plt[[1]]$layers[[2]]$geom_params$seed <- 42

    ggsave(
      filename = filename,
      plot = .plt,
      width = width,
      height = 4,
      units = "in",
      dpi = 600
    )
  }
)

# ------------------------------------------------------------------------------
# dot plot of RL lineage markers (including markers for control cells like
# endothelial cells, microglia, and oligodendrocytes); code taken from
# `rl_gene_expr.qmd` in thesis figures directory

# marker genes for the general cell types; selected from /.mounts/labs/pailab/private/icheong/CBL_scRNAseq/scrnaseq_Leo/utilities/cell_gene_mapping.csv
cell_type_markers <- c(
  "MKI67", # RL, GCP
  "WLS", # RL, UBC progenitors
  "EOMES", "LMX1A", "OTX2", # UBC
  "RBFOX3", "ATOH1", # GCP
  "PAX6", # differentiated cells (UBC/GC)
  "RELN", # GC
  "PDGFRA", "TNR", # OPC/oligodendrocytes
  "TREM2", # microglia
  "CLDN5", "ITM2A" # endothelial
)

.plt <- DotPlot(
  srat$rl,
  assay = "RNA",
  features = rev(cell_type_markers),
  group.by = "broad_annot"
) + 
  coord_flip() + 
  labs(x = "gene", y = "cell type") + 
  theme_classic2(base_size = 12) + 
  theme(
    axis.text.x = element_text(hjust = 1, angle = 30)
  )

ggsave(
  filename = "rl_lineage_markers_dotplot.png",
  plot = .plt,
  width = 6,
  height = 6,
  units = "in",
  dpi = 600
)

# ------------------------------------------------------------------------------
# scProportionTest for RL lineage; code taken from `rl_props.qmd` in thesis
# figures directory

# filter out control cells (endothelial cells, microglia, oligodendrocytes)
rl_prop_srat <- subset(
  x = srat$rl,
  subset = broad_annot %in% c("endothelial", "microglia", "oligodendrocyte/OPC"),
  invert = TRUE
)
rl_prop_srat$broad_annot <- fct_drop(rl_prop_srat$broad_annot)

# proportion of UBCs in human vs mouse RL
.plt <- cluster_barplot(
  object = rl_prop_srat,
  split.by = "broad_annot",
  group.by = "species",
  position = "fill"
) + 
  labs(fill = "cell type") + 
  scale_fill_manual(values = my_pals$rl_integ_annot) + 
  theme_classic2()
ggsave(
  filename = "rl_prop_bar.png",
  plot = .plt,
  width = 3,
  height = 3.75,
  units = "in",
  dpi = 600
)

# load scProportionTest results
scpt_res <- read_csv("/.mounts/labs/pailab/private/icheong/CBL_scRNAseq/results/integrated/20240822/rl_broad_annot_results.csv")

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
  clusters = fct_reorder(clusters, obs_log2FD, .desc = FALSE)
)

# plot scProportionTest results
scpt_plt <- ggplot(
  data = scpt_res,
  mapping = aes(
    x = obs_log2FD,
    y = clusters,
    xmin = boot_CI_2.5,
    xmax = boot_CI_97.5,
    colour = significance
  )
) + 
  geom_pointrange(size = 0.1) + 
  geom_vline(
    xintercept = c(log2FD_threshold, -log2FD_threshold),
    lty = "dashed"
  ) + 
  geom_vline(xintercept = 0, lty = "solid") + 
  labs(x = "observed log2(FD)", y = "cell type", colour = "enriched") + 
  scale_colour_manual(values = c(my_pals$species[1], "black", my_pals$species[2])) + 
  theme_classic2() + 
  theme(legend.position = "bottom")
ggsave(
  filename = "rl_permutation_plot.png",
  plot = scpt_plt,
  width = 4,
  height = 2.5,
  units = "in",
  dpi = 900
)

# ------------------------------------------------------------------------------
# pathways enriched in UBC_1 and UBC_2; code taken from `ubc_genes_regulons.qmd`
# in thesis figures directory

# load pathway enrichment results for all UBC subclusters
ubc_pathways <- map(
  .x = c("UBC_1", "UBC_2"),
  .f = \(clust) {
    ubc_pathways <- read_csv(file.path(
      "/.mounts/labs/pailab/private/llau/results/integrated/20240715",
      clust,
      "upregulated.csv"
    )) %>%
      # clean up the pathway name
      mutate(
        term_id = str_remove(string = term_id, pattern = "%.*")
      )

    return(ubc_pathways)
  }
) %>%
  setNames(nm = c("UBC_1", "UBC_2"))

# enriched pathways for all UBC subclusters
walk2(
  .x = ubc_pathways,
  .y = names(ubc_pathways),
  .f = \(pway, clust) {
    # total number of enriched pathways
    caption <- paste0(
      "(total # enriched pathways = ",
      nrow(filter(pway, p_value < 0.05)),
      ")"
    )

    # take top 10 pathways
    pway <- slice_min(
      pway,
      order_by = p_value,
      n = 10
    )

    # plot enriched pathways
    .plt <- ggplot(
      pway,
      aes(
        x = -log10(p_value),
        y = fct_reorder(.f = term_id, .x = p_value, .desc = TRUE)
      )
    ) + 
      geom_col(fill = my_pals$ubc_integ_clust[[clust]], width = 0.75) + 
      labs(
        caption = caption,
        x = "-log10(FDR)",
        y = "pathway term"
      ) + 
      scale_x_continuous(expand = expansion(mult = c(0, 0.05))) + 
      scale_y_discrete(label = scales::label_wrap(25)) + 
      theme_classic2() + 
      theme(plot.caption = element_text(size = rel(1)))

    ggsave(
      filename = paste0(clust, "_pathways.png"),
      plot = .plt,
      width = 4.5,
      height = 5,
      units = "in",
      dpi = 600
    )
  }
)

# ------------------------------------------------------------------------------
# regulon specificity score for UBC_1; code taken from `ubc_genes_regulons.qmd`
# in thesis figures directory

# load regulon specificity scores
rss <- read.csv(
  file = "/.mounts/labs/pailab/private/icheong/CBL_scRNAseq/results/pyscenic/20240704/cell_type_annot/aldinger_sepp_RL.rss.csv",
  row.names = 1
)
ubc_subclusters <- paste0("UBC_", c(0:5))

# rank top regulons for all UBCs
ubc_top_rss <- map(
  .x = ubc_subclusters,
  .f = \(cluster) {
    df <- data.frame(
      gene = colnames(rss),
      score = rss[cluster, ] %>% as.numeric(),
      rank = rank(-rss[cluster, ], ties.method = "random") %>% as.integer()
    ) %>%
      mutate(gene_label = case_when(
        rank <= 8 ~ gene,
        .default = ""
      )) %>%
      # change order for plotting
      arrange(desc(rank))

    return(df)
  }
) %>%
  setNames(nm = ubc_subclusters)

# RSS plots for UBC_1 and UBC_2 only
walk2(
  .x = ubc_top_rss[2:3],
  .y = names(ubc_top_rss[2:3]),
  .f = \(top_rss, clust) {
    .plt <- ggplot(
      top_rss,
      aes(x = rank, y = score)
    ) + 
      geom_point(colour = ifelse(
        top_rss$gene_label == "",
        "black",
        my_pals$ubc_integ_clust[[clust]]
      )) + 
      ggrepel::geom_text_repel(
        aes(label = gene_label),
        colour = my_pals$ubc_integ_clust[[clust]],
        size = 4,
        segment.colour = alpha(my_pals$ubc_integ_clust[[clust]], 0.5),
        min.segment.length = 0.1,
        max.overlaps = Inf,
        nudge_y = -0.005,
        show.legend = FALSE,
        seed = 100
      ) + 
      labs(x = "rank", y = "regulon specificity score") + 
      theme_classic2(base_size = 14)

    ggsave(
      filename = paste0(clust, "_rss_plot.png"),
      plot = .plt,
      width = 3,
      height = 4,
      units = "in",
      dpi = 600
    )
  }
)

# ------------------------------------------------------------------------------
# SOX4 and SOX11 regulon expression in the RL lineage

# load the human-only Seurat object (Aldinger/Sepp)
srat$human_rl <- readRDS("/.mounts/labs/pailab/private/icheong/CBL_scRNAseq/results/pyscenic/20240704/integ_human_srat.rds")

# import regulon expression values
human_rl_regulon_expr <- read.csv(
  "/.mounts/labs/pailab/private/icheong/CBL_scRNAseq/results/pyscenic/20240926/aldinger_sepp_RL_regulon_expr.csv",
  row.names = 1,
  check.names = FALSE
)

# add regulon expression to metadata
human_rl_regulon_expr <- human_rl_regulon_expr[rownames(srat$human_rl[[]]), ]
srat$human_rl@meta.data <- human_rl_regulon_expr

# set idents for dot plot
Idents(srat$human_rl) <- "cell_type_annot"

# plot SOX4 and SOX11 regulon expression
walk(
  .x = c("SOX4", "SOX11"),
  .f = \(regulon) {
    .plt <- FeaturePlot(
      srat$human_rl,
      features = paste0("Regulon(", regulon, ")"),
      cols = scales::pal_viridis()(100),
      order = TRUE
    )

    ggsave(
      filename = paste0("human_rl_", regulon, ".png"),
      plot = .plt,
      width = 4,
      height = 4,
      units = "in",
      dpi = 600
    )
  }
)

# dot plot of regulon expression
regs_to_plot <- c("OTX2", "SOX4", "BACH1", "SOX11")
.plt <- DotPlot(
  srat$human_rl,
  features = rev(paste0("Regulon(", regs_to_plot, ")")),
  idents = paste0("UBC_", c(0:5)),
  group.by = "cell_type_annot",
  scale = TRUE
) + 
  labs(x = "regulons", y = "cell type") + 
  coord_flip() + 
  theme_classic2()
ggsave(
  filename = "human_rl_regulon_dotplot.png",
  plot = .plt,
  width = 6,
  height = 3,
  units = "in",
  dpi = 600
)

# ------------------------------------------------------------------------------
# medulloblastoma UMAPs

walk(
  .x = c("seurat_clusters", "subtype"),
  .f = \(group_by) {
    .plt <- DimPlot(
      srat$mb,
      reduction = "umap",
      group.by = group_by,
      label = TRUE,
      repel = TRUE
    ) + 
      labs(title = NULL) + 
      theme_classic2() + 
      NoLegend()

    if (group_by == "seurat_clusters") {
      filename <- "mb_clusters.png"
    } else if (group_by == "subtype") {
      .plt <- .plt + 
        scale_colour_manual(values = my_pals$mb_subtype)
      filename <- "mb_subtypes.png"
    }

    ggsave(
      filename = filename,
      plot = .plt,
      width = 4,
      height = 4,
      units = "in",
      dpi = 600
    )
  }
)

# ------------------------------------------------------------------------------
# medulloblastoma volcano plot; code taken from 2024-09-17 notes

# load differential gene expression results
mb_diff_exp <- read_csv("/.mounts/labs/pailab/private/icheong/CBL_scRNAseq/results/tumour/Vladoiu/20240917/mb_subtype_diff_exp.csv")

# set genes with `p_val_adj` of 0 to 1e-300 (the minimum non-zero `p_val_adj` is
# around 3e-294 right now)
min_adj_p <- filter(mb_diff_exp, p_val_adj != 0) %>% pull(p_val_adj) %>% min()
mb_diff_exp <- mutate(
  mb_diff_exp,
  p_val_adj = case_when(
    p_val_adj < min_adj_p ~ 1e-300,
    .default = p_val_adj
  )
)

# add columns for volcano plot
mb_diff_exp_tfs <- map(
  .x = mb_subtypes,
  .f = \(subtype) {
    de <- filter(mb_diff_exp, cluster == subtype & gene %in% tfs) %>%
      mutate(
        direction = case_when(
          avg_log2FC > 0 & p_val_adj < 0.05 ~ "up",
          avg_log2FC < 0 & p_val_adj < 0.05 ~ "down",
          .default = "n.s."
        )
      ) %>%
      arrange(desc(avg_log2FC), p_val_adj)

    # label top 10 upreg and top 10 downreg TFs
    genes_to_label <- c(
      head(de$gene, n = 10),
      tail(de$gene, n = 10),
      "SOX4",
      "SOX11",
      "OTX2",
      "EOMES",
      "LMX1A"
    )
    de <- mutate(
      de,
      gene_label = case_when(
        gene %in% genes_to_label ~ gene,
        .default = ""
      ),
    )

    return(de)
  }
) %>%
  setNames(nm = mb_subtypes)

# volcano plot for differentially expressed TFs in G4 MB
.plt <- make_volcano(
  data = mb_diff_exp_tfs$G4,
  log_fc = avg_log2FC,
  log_pval = -log10(p_val_adj),
  direction = direction,
  gene_labels = gene_label,
  log_fc_thresh = 0,
  log_pval_thresh = -log10(0.05),
  seed = 4
) + 
  labs(
    x = "log2(FC)",
    y = "-log10(adj. p-value)",
    colour = "Group 4\nTFs"
  ) + 
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.10))) + 
  scale_colour_manual(
    values = c("down" = "blue", "n.s." = "black", "up" = "red"),
    guide = guide_legend(override.aes = list(size = 2))
  ) + 
  theme_classic2() + 
  theme(plot.background = element_rect(fill = "white", colour = NA))
# set text label size
.plt$layers[[2]]$aes_params$size <- 2.5

ggsave(
  filename = "mb_G4_TF_volcano.png",
  plot = .plt,
  width = 5,
  height = 4,
  units = "in",
  dpi = 600
)

