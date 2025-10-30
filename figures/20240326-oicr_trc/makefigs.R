# ==============================================================================
# These figures were generated for the OICR Translational Research Conference in
# March 2024. (Modified from the November 2023 CEEHRC conference.)
# ==============================================================================

# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(tidyverse)
library(patchwork)
library(ggrepel)
library(Seurat)
# library(SeuratDisk)

# >>> load Seurat

# Seurat objects for import
srat_names <- c("aldinger", "vladoiu_mb")
srat_rds <- c(
  # Aldinger human RL from Liam
  "/isilon/CBL_scRNAseq-archived/data/human/Aldinger/glutamatergic_dev_Liam.RDS",
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

# <<< Seurat

# >>> load pySCENIC results

# add regulon expression to Seurat metadata
regulon_expr <- map(
  .x = c(
    "../20231113-ceehrc_conference/aldinger_rl_regulon-20230725.csv",
    "../20231113-ceehrc_conference/vladoiu_mb_regulon-20231018.csv"
  ),
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
        gene_label = case_when(gene %in% genes ~ gene, TRUE ~ "")
      ) %>% 
      arrange(gene_label) # order so the labelled genes get plotted on top
    
    .plt <- ggplot(df, aes(x = rank, y = rss, label = gene_label)) + 
      geom_point(color = ifelse(df$gene_label == "", "black", "#E41A1C"), size = 1) + 
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
