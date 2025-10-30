# ==============================================================================
# Supplementary figure of manuscript showing the top regulons in each UBC
# cluster. This script should be run using the scrnaseq_env conda environment.
# ==============================================================================

library(tidyverse)
library(patchwork)

source("./utils.R")
source("../../software/utilities/pyscenic_regulons.R")

out_dir <- file.path(
  "/.mounts/labs/pailab/private/projects/HumanMouseUBC/figures/SuppFig_ubc_regulons",
  format(Sys.Date(), "%Y%m%d")
)
if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

# colour palette
my_pals <- get_custom_pals()

# ------------------------------------------------------------------------------
# RSS plots

# load regulon specificity scores
rss <- read.csv(
  file = "/.mounts/labs/pailab/private/icheong/CBL_scRNAseq/results/pyscenic/20240704/cell_type_annot/aldinger_sepp_RL.rss.csv",
  row.names = 1
)

# rank top regulons for all UBCs
ubc_top_rss <- map(
  .x = paste0("UBC_", c(0:5)),
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
  setNames(nm = paste0("iUBC", c(0:5)))

# RSS plots
rss_plots <- map2(
  .x = ubc_top_rss,
  .y = names(ubc_top_rss),
  .f = \(top_rss, clust) {
    .plt <- ggplot(
      top_rss,
      aes(x = rank, y = score)
    ) + 
      geom_point(colour = ifelse(
        top_rss$gene_label == "",
        "black",
        my_pals$ubc_integ_clust2[[clust]]
      )) + 
      ggrepel::geom_text_repel(
        aes(label = gene_label),
        colour = my_pals$ubc_integ_clust2[[clust]],
        size = 2.5,
        segment.colour = alpha(my_pals$ubc_integ_clust2[[clust]], 0.5),
        min.segment.length = 0.1,
        max.overlaps = Inf,
        nudge_x = 5,
        show.legend = FALSE,
        seed = 100
      ) + 
      labs(title = clust, x = "rank", y = "regulon specificity score") + 
      theme_classic2()

    ggsave(
      filename = paste0(clust, "_rss_plot.pdf"),
      plot = .plt,
      path = out_dir,
      width = 3,
      height = 4,
      units = "in"
    )

    return(.plt)
  }
)

# ------------------------------------------------------------------------------
# final figure

.plt <- wrap_plots(rss_plots, nrow = 2) 

walk(
  .x = c("png", "pdf"),
  .f = \(device) {
    ggsave(
      filename = paste0("combined.", device),
      plot = .plt,
      path = out_dir,
      width = 9,
      height = 8,
      units = "in",
      dpi = 600
    )
  }
)
