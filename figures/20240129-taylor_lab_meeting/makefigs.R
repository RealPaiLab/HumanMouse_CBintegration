# ==============================================================================
# These figures were generated for the Taylor Lab meeting presentation in
# January 2024.
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

# RSS
all_rss <- map(
  .x = list(
    rl_rss = "/CBL_scRNAseq/results/pyscenic/20230725/integ/aldinger_RL.rss.csv",
    mb_rss = "/CBL_scRNAseq/results/pyscenic/20231018/vladoiu_mb.rss.csv"
  ),
  .f = \(x) {read.csv(x, row.names = 1)}
)

# <<< pySCENIC

# colour palette for species
species_cols <- c(hcl(h = 15, c = 100, l = 65), "grey")

# import functions
source("/CBL_scRNAseq/software/utilities/plotting.R")

# ------------------------------------------------------------------------------
# RL-SVZ bar plot

.plt <- cluster_barplot(
  all_srat$cca_RL,
  split.by = "species",
  width = 0.6,
  filter_data = "seurat_clusters %in% c(8)"
) + 
  labs(x = "clusters", fill = NULL) + 
  scale_fill_manual(labels = c("mouse", "human"), values = rev(species_cols)) + 
  theme_classic() + 
  theme(
    legend.position = "right",
    axis.text = element_text(colour = "black", size = 11),
    axis.ticks = element_line(colour = "black")
  )

ggsave(
  "rlsvz_bar.png",
  plot = .plt,
  width = 2.5,
  height = 4,
  units = "in",
  dpi = 600
)

# ------------------------------------------------------------------------------
# RL-SVZ/tumour RSS plot

cells_to_plot <- tribble(
  ~rss, ~cluster, ~filename,
  "rl_rss", "7-Homol UBC", "rss_7homol.png",
  "rl_rss", "19-NonHomol UBC", "rss_19nonhomol.png",
  "rl_rss", "20-NonHomol UBC", "rss_20nonhomol.png",
  "rl_rss", "8-Human GCP/RL_SVZ", "rss_8rlsvz.png",
  "mb_rss", "6", "rss_mb6.png",
  "mb_rss", "8", "rss_mb8.png",
)

pwalk(
  .l = cells_to_plot,
  .f = \(rss, cluster, filename) {
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
        gene_label = case_when(rank <= 8 ~ gene, TRUE ~ NA_character_)
      )
    
    .plt <- ggplot(df, aes(x = rank, y = rss, label = gene_label)) + 
      geom_point(color = ifelse(is.na(df$gene_label), "black", "#E41A1C"), size = 1) + 
      geom_text_repel(color = "#E41A1C", fontface = "bold") + 
      labs(x = "Rank", y = "Regulon Specificity Score (RSS)", title = cluster) + 
      theme_classic() + 
      theme(
        axis.text = element_text(size = 8, colour = "black"),
        axis.ticks = element_line(colour = "black"),
        axis.title = element_text(size = 10)
      )
    
    ggsave(
      filename = filename,
      plot = .plt,
      width = 2.5,
      height = 4,
      units = "in",
      dpi = 600
    )
  }
)

# ------------------------------------------------------------------------------
# BRCA1 expression in developing hindbrain and tumour

.plt <- FeaturePlot(
  all_srat$aldinger,
  features = "BRCA1",
  order = TRUE
)
ggsave(
  filename = "brca1_expr_dev.png",
  plot = .plt,
  width = 5,
  height = 4,
  units = "in",
  dpi = 600
)

.plt <- FeaturePlot(
  all_srat$vladoiu_mb,
  features = "BRCA1",
  order = TRUE
)
ggsave(
  filename = "brca1_expr_mb.png",
  plot = .plt,
  width = 5,
  height = 4,
  units = "in",
  dpi = 600
)
