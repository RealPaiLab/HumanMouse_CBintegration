# ==============================================================================
# These figures were generated for a grant proposal, requested by Shraddha on
# 2024-03-19.
# ==============================================================================

# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(tidyverse)
library(Seurat)
library(ggrepel)

# ------------------------------------------------------------------------------
# UMAP of Aldinger data with Hendrikse annotations

hend_srat <- readRDS("/isilon/CBL_scRNAseq-archived/data/human/Aldinger/glutamatergic_dev_Liam.RDS")

.plt <- DimPlot(hend_srat, reduction = "umap", group.by = "new_cell_type", label = FALSE) + 
  labs(title = "Human rhombic lip-derived cells") + 
  theme(
    axis.ticks = element_blank(),
    axis.text = element_blank()
  )
ggsave(
  filename = "hendrikse_umap.png",
  plot = .plt,
  width = 5,
  height = 4,
  units = "in",
  dpi = 600
)

# ------------------------------------------------------------------------------
# pySCENIC RSS score

all_rss <- map(
  .x = list(
    rl_rss = "/CBL_scRNAseq/results/pyscenic/20230725/hendrikse/aldinger_RL.rss.csv",
    mb_rss = "/CBL_scRNAseq/results/pyscenic/20231018/vladoiu_mb.rss.csv"
  ),
  .f = \(x) {read.csv(x, row.names = 1)}
)

# which RSS clusters to plot
cells_to_plot <- tribble(
  ~rss, ~cluster, ~filename,
  "rl_rss", "RL-SVZ", "rlsvz_rss.png",
  "mb_rss", "6", "mb6_rss.png",
  "mb_rss", "8", "mb8_rss.png",
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
        gene_label = case_when(rank <= 8 ~ gene, TRUE ~ "")
      )
    
    .plt <- ggplot(df, aes(x = rank, y = rss, label = gene_label)) + 
      geom_point(color = ifelse(df$gene_label == "", "black", "#E41A1C"), size = 0.5) + 
      geom_text_repel(
        color = "#E41A1C",
        fontface = "bold",
        segment.size = 0.4,
        bg.colour = "white",
        bg.r = 0.025
      ) + 
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
      width = 3,
      height = 3,
      units = "in",
      dpi = 600
    )
  }
)
