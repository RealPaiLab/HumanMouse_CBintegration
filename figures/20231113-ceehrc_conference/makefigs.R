# ==============================================================================
# These figures were generated for the CEEHRC Canadian Conference on Epigenetics
# poster in November 2023.
# ==============================================================================

# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(tidyverse)
library(patchwork)
library(ggrepel)
library(Seurat)
library(SeuratDisk)

# >>> load Seurat

# Seurat objects for import
srat_names <- c("aldinger", "aldinger_full", "vladoiu", "cca_RL", "cca_full", "vladoiu_mb")
srat_rds <- c(
  # Aldinger human RL from Liam
  "/isilon/CBL_scRNAseq-archived/data/human/Aldinger/glutamatergic_dev_Liam.RDS",
  # full Aldinger dataset
  "/isilon/CBL_scRNAseq-archived/data/human/Aldinger/seurat.rds",
  # Vladoiu mouse full cerebellum
  "/CBL_scRNAseq/results/mouse/Vladoiu/merged_seurat.rds",
  # CCA RL only
  "/CBL_scRNAseq/results/integrated/vladoiu_liam_RL.rds",
  # CCA full cerebellum
  "/CBL_scRNAseq/results/integrated/20230126/without_future/aldinger_vladoiu_cca.rds",
  # medulloblastoma
  "/CBL_scRNAseq/results/tumour/Vladoiu/20230510/mb_mnn.rds"
)

# import Seurat objects
all_srat <- map(
  .x = srat_rds,
  .f = readRDS
)
names(all_srat) <- srat_names
print(all_srat)

# convert to factor
all_srat$cca_RL$species <- fct_relevel(all_srat$cca_RL$species, "mouse", "human")

# <<< Seurat

# >>> load pySCENIC results

# add regulon expression to Seurat metadata
regulon_expr <- map(
  .x = c("aldinger_rl_regulon-20230725.csv", "vladoiu_mb_regulon-20231018.csv"),
  .f = \(x) {
    expr <- read.csv(x, row.names = 1)
    expr <- rename_with(expr, ~ str_remove(.x, ".$"), starts_with("Regulon"))
    return(expr)
  }
)
names(regulon_expr) <- c("rl", "mb")

all_srat$aldinger@meta.data <- regulon_expr$rl
all_srat$vladoiu_mb@meta.data <- regulon_expr$mb

# RSS
all_rss <- map(
  .x = list(
    rl_rss = "/CBL_scRNAseq/results/pyscenic/20230725/integ/aldinger_RL.rss.csv",
    mb_rss = "/CBL_scRNAseq/results/pyscenic/20231018/vladoiu_mb.rss.csv"
  ),
  .f = \(x) {read.csv(x, row.names = 1)}
)

# <<< pySCENIC

# import functions
# source("/CBL_scRNAseq/software/utilities/cluster_barplot.R")
source("/CBL_scRNAseq/software/utilities/plotting.R")

# colour palette for species
species_cols <- c("#0099AD", "#FFC857")
# species_cols <- c("#0099AD", "#FED9B7")

# size for axis.title
umap_ax_size = 11

# ------------------------------------------------------------------------------
# human/mouse UMAP pre-integration

human_cells <- rownames(all_srat$cca_full@meta.data)[all_srat$cca_full$species == "human"]

# UMAP of full Aldinger dataset
human_full.plt <- DimPlot(
  all_srat$aldinger_full,
  reduction = "umap",
  cells.highlight = rownames(all_srat$aldinger_full[[]]),
  sizes.highlight = 0.01,
  label = FALSE
) + 
  NoLegend() + 
  scale_colour_manual(labels = "human", values = species_cols[1]) +
  labs(title = "Human cerebellum")
  # theme(axis.text = element_blank(), axis.ticks = element_blank())

# UMAP of mouse RL-derived cells
mouse_full.plt <- DimPlot(
  all_srat$vladoiu,
  reduction = "umap",
  cells.highlight = rownames(all_srat$vladoiu[[]]),
  sizes.highlight = 0.01,
  label = FALSE,
) + 
  NoLegend() + 
  scale_colour_manual(labels = "mouse", values = species_cols[2]) + 
  labs(title = "Mouse cerebellum")
  # theme(axis.text = element_blank(), axis.ticks = element_blank())

# UMAP of human + mouse integration entire cerebellum
integ_full.plt <- DimPlot(
  all_srat$cca_full,
  reduction = "umap",
  group.by = "species",
  pt.size = 0.01,
  cells.highlight = human_cells,
  sizes.highlight = 0.01,
  label = FALSE,
  raster = FALSE
) + 
  NoLegend() + 
  scale_colour_manual(labels = c("mouse", "human"), values = rev(species_cols)) + 
  labs(title = "Human + mouse integration")
  # theme(axis.text = element_blank(), axis.ticks = element_blank())

# UMAP of human + mouse integration RL only
integ_rl.plt <- DimPlot(
  all_srat$cca_RL,
  reduction = "umap",
  group.by = "species",
  pt.size = 0.01,
  cells.highlight = human_cells,
  sizes.highlight = 0.01,
  label = FALSE
) + 
  NoLegend() + 
  scale_colour_manual(labels = c("mouse", "human"), values = rev(species_cols)) + 
  labs(title = "Human + mouse rhombic lip")
  # theme(axis.text = element_blank(), axis.ticks = element_blank())

layout <- c(
  area(t = 1, l = 1, r = 7),
  area(t = 2, l = 1, r = 7),
  area(t = 1, l = 10, r = 16),
  area(t = 2, l = 10, r = 16)
)
.plt <- human_full.plt + mouse_full.plt + integ_full.plt + integ_rl.plt + 
  plot_layout(byrow = FALSE, design = layout) + 
  plot_annotation(tag_levels = "A") & 
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 14),
    plot.tag = element_text(vjust = 0.5),
    plot.tag.position = c(0, 0.975),
    axis.text = element_blank(),
    axis.ticks = element_blank()
  )
ggsave(
  "integration.png",
  plot = .plt,
  width = 8,
  height = 6,
  units = "in",
  dpi = 1200
)

# ------------------------------------------------------------------------------
# tumour dataset

mb_umap.plt <- DimPlot(
  all_srat$vladoiu_mb,
  group.by = "subtype",
  cols = RColorBrewer::brewer.pal(3, "Dark2"),
  label = TRUE,
  repel = TRUE
) + 
  NoLegend() + 
  labs(title = "MB tumour") + 
  theme(axis.text = element_blank(), axis.ticks = element_blank())
ggsave(
  "mb_subtype_umap.png",
  plot = mb_umap.plt,
  width = 3,
  height = 3,
  units = "in",
  dpi = 1200
)

# ------------------------------------------------------------------------------
# human-specific UBCs

# UBC cell IDs
ubcs <- purrr::map(
  .x = c(7, 19, 20),
  .f = \(x) {
    rownames(all_srat$cca_RL@meta.data)[all_srat$cca_RL$seurat_clusters %in% x]
  }
)

# EOMES expression
eomes.plt <- FeaturePlot(
  all_srat$cca_RL,
  features = "EOMES",
  min.cutoff = "q10",
  max.cutoff = "q90",
  order = TRUE
) + 
  labs(title = "*EOMES* (UBC marker)") + 
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    plot.title = ggtext::element_markdown()
  )

# highlight UBC clusters
ubcs.plt <- DimPlot(
  all_srat$cca_RL,
  reduction = "umap",
  group.by = "seurat_clusters",
  cells.highlight = ubcs,
  cols.highlight = rev(RColorBrewer::brewer.pal(n = 3, name = "Dark2")),
  sizes.highlight = 0.01,
  label = TRUE,
  repel = TRUE
) + 
  NoLegend() + 
  labs(title = "Three UBC clusters") + 
  theme(axis.text = element_blank(), axis.ticks = element_blank())

# bar plot - species per cluster
ubc_species.plt <- cluster_barplot(
  all_srat$cca_RL,
  split.by = "species",
  width = 0.6,
  filter_data = "seurat_clusters %in% c(7, 19, 20)"
) + 
  labs(x = "UBC clusters", y = "Number of cells", fill = NULL) + 
  scale_fill_manual(labels = c("mouse", "human"), values = rev(species_cols)) + 
  theme_classic() + 
  theme(
    legend.position = c(0.75, 0.8),
    axis.text = element_text(colour = "black"),
    axis.ticks = element_line(colour = "black")
  )

.plt <- eomes.plt + ubcs.plt + wrap_elements(full = ubc_species.plt) + 
  plot_layout(ncol = 3) + 
  plot_annotation(tag_levels = "A") & 
  theme(
    axis.title = element_text(size = 11),
    plot.tag = element_text(face = "bold", size = 14, vjust = 1),
    plot.tag.position = c(0, 1)
  )
ggsave(
  "human_specific_ubcs.png",
  plot = .plt,
  width = 9,
  height = 3,
  units = "in",
  dpi = 1200
)

# ------------------------------------------------------------------------------
# differential expression, enriched TFs

# import differential gene expression results
de_genes <- read_csv("/CBL_scRNAseq/results/human/Aldinger/20230205/de_genes.csv")

# set threshold for differential expression
logfc_threshold <- 1

# add labels
# human-specific UBC genes get swapped for the plot
de_genes <- mutate(
  de_genes,
  avg_log2FC = -1 * avg_log2FC,
  diff_exp = case_when(
    (avg_log2FC < 0 & p_val_adj < 0.05) ~ "down",
    (avg_log2FC > 0 & p_val_adj < 0.05) ~ "up",
    TRUE ~ "ns"
  ) %>% factor(levels = c("up", "down", "ns")),
  gene_label = case_when(abs(avg_log2FC) > 1 ~ gene)
)


# make volcano plot
volcano.plt <- make_volcano(
  data = de_genes,
  log_fc = avg_log2FC,
  log_pval = -log10(p_val_adj),
  direction = diff_exp,
  label_num_genes = FALSE
) + 
  geom_text_repel(aes(label = gene_label), size = 3) + 
  geom_hline(
    yintercept = -log10(0.05),
    colour = "grey50",
    linewidth = 0.5,
    linetype = "dashed"
  ) + 
  labs(
    x = expression("log"[2]*"(fold change)"),
    y = expression("-log"[10]*"(adjusted p-value)")
  ) + 
  scale_colour_manual(
    values = c(RColorBrewer::brewer.pal(3, "Set1")[1:2], "black")
  ) + 
  theme_classic() + 
  theme(
    axis.text = element_text(colour = "black"),
    axis.ticks = element_line(colour = "black"),
    legend.position = "none"
  )
ggsave(
  "ubc_diff_exp.png",
  plot = volcano.plt,
  width = 3,
  height = 3,
  units = "in",
  dpi = 1200
)



# import list of AME TFs
enr_motifs <- read_csv("/CBL_scRNAseq/results/human/Aldinger/20230224/enriched_motifs_nonhomol.csv")

# subset GFI1, GFI1B, PRDM6, MYCN, ZIC1
enr_motifs_sub <- enr_motifs %>% 
  mutate(
    motif_id = str_remove_all(string = motif_id, pattern = "_.*")
  ) %>% 
  filter(
    motif_id %in% c("GFI1", "GFI1B", "PRDM6", "MYCN", "ZIC1")
  ) %>% 
  select(motif_id, adj.pvalue)

tf_motifs.plt <- memes::plot_ame_heatmap(enr_motifs_sub) + 
  # add numbers
  geom_text(
    aes(x = motif_id, y = "All Regions", label = round(-log10(adj.pvalue), 1)),
    data = enr_motifs_sub
  ) + 
  labs(
    x = "Enrichment of known MB drivers",
    tag = "C",
    fill = "-log<sub>10</sub>(adj.<br>p-value)"
  ) + 
  scale_x_discrete(expand = expansion(mult = 0), limits = rev) + 
  scale_y_discrete(expand = expansion(mult = 0)) + 
  coord_flip() + 
  theme(
    axis.ticks.x = element_blank(),
    axis.text.x = element_blank(),
    legend.title = ggtext::element_markdown(),
    plot.tag = element_text(face = "bold", size = 14)
  )
ggsave(
  "key_enr_motifs.png",
  plot = tf_motifs.plt,
  width = 2.75,
  height = 3.5,
  units = "in",
  dpi = 1200
)

# ------------------------------------------------------------------------------
# UBC pySCENIC plots

rss.plt <- pmap(
  .l = list(
    list("7-Homol UBC", "19-NonHomol UBC", "20-NonHomol UBC"),
    list("7-Common UBC", "19-Human specific UBC", "20-Human specific UBC")
  ),
  .f = \(cluster, name, rss = all_rss$rl_rss) {
    rss <- rss[cluster, ]
    
    df <- data.frame(
      gene = names(rss[cluster, ]),
      # extract RSS score as vector
      rss = rss[cluster, ] %>% as.numeric(),
      # get rank, negative sign needed to rank the highest RSS as 1
      rank = rank(-rss[cluster, ], ties.method = "random") %>% as.integer()
    ) %>% 
      mutate(
        gene_label = case_when(rank <=5 ~ gene, TRUE ~ NA_character_)
      )
    
    .plt <- ggplot(df, aes(x = rank, y = rss, label = gene_label)) + 
      geom_point(color = ifelse(is.na(df$gene_label), "black", "#377EB8"), size = 1) + 
      geom_text_repel(nudge_x = 10, color = "#377EB8", size = 3) + 
      labs(x = "Rank", y = "Regulon Specificity Score (RSS)") + 
      theme_classic() + 
      theme(
        axis.text = element_text(size = 8, colour = "black"),
        axis.ticks = element_line(colour = "black"),
        axis.title = element_text(size = 10),
        plot.title = element_text(hjust = 0.5, face = "bold")
      )
    
    return(.plt)
  }
)
.plt <- wrap_plots(rss.plt, ncol = 3)
ggsave(
  "ubc_rss.png",
  plot = .plt,
  width = 4.5,
  height = 3,
  units = "in",
  dpi = 1200
)

# ------------------------------------------------------------------------------
# BRCA1 pySCENIC plots

# plot pySCENIC RSS
genes <- c("BRCA1")
rss.plt <- pmap(
  .l = list(
    list("rl_rss", "mb_rss", "mb_rss"),
    list("8-Human GCP/RL_SVZ", "6", "8")
  ),
  .f = \(rss, cluster) {
    rss <- all_rss[[rss]]
    rss <- rss[rownames(rss) %in% cluster, ]
    # rownames(rss) <- str_replace_all(rownames(rss), "[:blank:]|/", "_")
    
    df <- data.frame(
      gene = names(rss[cluster, ]),
      # extract RSS score as vector
      rss = rss[cluster, ] %>% as.numeric(),
      # get rank, negative sign needed to rank the highest RSS as 1
      rank = rank(-rss[cluster, ], ties.method = "random") %>% as.integer()
    ) %>% 
      mutate(
        gene_label = case_when(gene %in% genes ~ gene, TRUE ~ NA_character_)
      )
    
    .plt <- ggplot(df, aes(x = rank, y = rss, label = gene_label)) + 
      geom_point(color = ifelse(is.na(df$gene_label), "black", "#E41A1C"), size = 1) + 
      geom_text_repel(color = "#E41A1C", fontface = "bold") + 
      labs(x = "Rank", y = "Regulon Specificity\nScore (RSS)") + 
      theme_classic() + 
      theme(
        axis.text = element_text(size = 8, colour = "black"),
        axis.ticks = element_line(colour = "black"),
        axis.title = element_text(size = 10)
      )
    
    return(.plt)
  }
)

# plot BRCA regulon expression
brca1_expr.plt <- pmap(
  .l = list(
    list(all_srat$aldinger, all_srat$vladoiu_mb),
    list("rl", "mb")
  ),
  .f = \(srat, sample) {
    .plt <- FeaturePlot(
      object = srat,
      features = "Regulon.BRCA1",
      order = TRUE,
      min.cutoff = "q10"
    ) + 
      labs(title = NULL) + 
      scale_colour_viridis_c(breaks = c(0.1, 0.2)) + 
      theme(
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_text(size = umap_ax_size)
      )
    
    return(.plt)
  }
)

# plot labelled UMAPs
rl_labs.plt <- DimPlot(
  all_srat$aldinger,
  group.by = "integ_clusters",
  label = TRUE,
  label.size = 3,
  repel = TRUE
) + 
  NoLegend() +
  labs(title = NULL) + 
  scale_colour_manual(values = unname(pals::trubetskoy())) + 
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_text(size = umap_ax_size)
  )
mb_labs.plt <- DimPlot(
  all_srat$vladoiu_mb,
  group.by = "seurat_clusters",
  label = TRUE,
  label.size = 3,
  repel = TRUE
) + 
  NoLegend() + 
  labs(title = NULL) + 
  scale_colour_manual(values = unname(pals::polychrome())) + 
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_text(size = umap_ax_size)
  )

# save plots
# design <-
#   "AAACFF
#    BBDEGG"
design <-
  "AAAABBBB
   AAAABBBB
   AAAABBBB
   #CC#DDEE
   #CC#DDEE
   FFFFGGGG
   FFFFGGGG
   FFFFGGGG"
.plt <- 
  wrap_elements(full = rl_labs.plt) + #A
  wrap_elements(full = mb_labs.plt) + #B
  rss.plt[[1]] + #C
  rss.plt[[2]] + #D
  rss.plt[[3]] + #E
  brca1_expr.plt[[1]] + #F
  brca1_expr.plt[[2]] + #G
  plot_annotation(tag_levels = "A") + 
  plot_layout(design = design) & 
  theme(plot.tag = element_text(face = "bold", size = 12))
ggsave(
  "brca1.png",
  plot = .plt,
  height = 9,
  width = 8,
  units = "in",
  dpi = 1200
)

# ------------------------------------------------------------------------------
# CRISPR target prioritization plot (DepMap)

source("prioritize_crispr_candidates.R") %>% suppressMessages()