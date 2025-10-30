# ==============================================================================
# These figures were generated for my committee meeting on 2024-02-26.
# ==============================================================================

# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(tidyverse)
library(patchwork)
library(ggrepel)

# load pathway enrichment results from human-specific/common UBCs
ubc_enr_csv <- c(
  "/CBL_scRNAseq/results/human/Aldinger/20240218/human_specific_ubc.csv",
  "/CBL_scRNAseq/results/human/Aldinger/20240218/common_ubc.csv"
)
ubc_enr <- map(
  .x = ubc_enr_csv,
  .f = read_csv
) %>% 
  `names<-`(c("human_specific", "common"))

# load list of mutations in BRCA1
brca1_mut <- read_csv("/CBL_scRNAseq/results/tumour/ICGC/20240225/BRCA1.csv") %>% 
  mutate(subgroup = factor(.$subgroup, levels = c("WNT", "SHH", "Group 3", "Group 4")))

# load BRCA1 CRISPRi differential expression data
brca1_kd <- data.table::fread(
  "/isilon/CBL_scRNAseq-archived/data/src/brca1_crispri/rnaseq/20240202_BRCA1-KD_KDvsNT_differential_expression_results_STAR_limma_trend.txt"
)

# load plotting functions
# source("/CBL_scRNAseq/software/utilities/plotting.R")

# ------------------------------------------------------------------------------
# plot top enriched pathways

pwalk(
  .l = list(
    enr_pathways = ubc_enr,
    cell_type = names(ubc_enr),
    fill_col = RColorBrewer::brewer.pal(3, "Set1")[1:2]
  ),
  .f = \(enr_pathways, cell_type, fill_col) {
    enr_pathways <- enr_pathways %>% 
      arrange(p_value) %>% 
      head(10) %>% 
      mutate(
        log_p_value = -log10(p_value)
      )
    
    .plt <- ggplot(
      enr_pathways,
      aes(x = fct_reorder(term_name, log_p_value),
          y = log_p_value)
    ) + 
      geom_col(fill = fill_col, width = 0.6) + 
      labs(x = "pathway", y = "-log10(adj. p-value)") + 
      scale_x_discrete(labels = scales::label_wrap(30)) + 
      scale_y_continuous(expand = expansion(mult = c(0, 0.05))) + 
      coord_flip() + 
      theme_classic() + 
      theme(
        axis.ticks = element_line(colour = "black"),
        axis.text = element_text(colour = "black")
      )
    
    ggsave(
      filename = paste0(cell_type, "_ubc_pathways.png"),
      plot = .plt,
      width = 3.5,
      height = 3.5,
      units = "in",
      dpi = 600
    )    
  }
)

# ------------------------------------------------------------------------------
# heatmap of shared regulons in development and MB

# see `<pyscenic software directory>/shared_regulons_dev_mb.R`

# ------------------------------------------------------------------------------
# mutations in BRCA1 (ICGC data)

# how many unique patients with BRCA1 mutation?
length(unique(brca1_mut$icgc_donor_id))

subtype_cols <- pals::brewer.set2(5)[2:5]

.plt <- ggplot(brca1_mut, aes(x = seq_type, fill = subgroup)) + 
  geom_bar(position = position_dodge()) + 
  labs(x = "gene region", y = "number of mutations") + 
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) + 
  scale_fill_manual(values = subtype_cols) + 
  theme_classic() + 
  theme(axis.ticks = element_line(colour = "black"),
        axis.text = element_text(colour = "black"))

ggsave(
  filename = "num_brca1_mutations.png",
  plot = .plt,
  width = 6,
  height = 4,
  units = "in",
  dpi = 600
)

# ------------------------------------------------------------------------------
# BRCA1 CRISPRi knockdown volcano plots

diff_gene_list <- c('PAX6', 'MKI67', 'RBFOX3', 'RELN', 'LMX1A', 'EOMES')
mb_reg_list <- read_csv("/CBL_scRNAseq/results/pyscenic/20231019/BRCA1.csv") %>% 
  pull(gene)

brca1_kd <- brca1_kd %>% 
  mutate(
    change = case_when(
      P.Value <= 0.05 & abs(logFC) >= 0.5 ~ "sig_dir",
      P.Value <= 0.05 & abs(logFC) < 0.5 ~ "sig_only",
      P.Value > 0.05 & abs(logFC) >= 0.5 ~ "dir_only",
      TRUE ~ "none"
    ) %>% factor(levels = c("sig_dir", "sig_only", "dir_only", "none")),
    sig_genes = case_when(
      SYMBOL %in% head(.$SYMBOL, 15) ~ SYMBOL,
      TRUE ~ NA_character_
    ),
    diff_genes = case_when(
      SYMBOL %in% diff_gene_list ~ SYMBOL,
      TRUE ~ NA_character_
    ),
    reg_genes = case_when(
      SYMBOL %in% mb_reg_list ~ SYMBOL,
      TRUE ~ NA_character_
    )
  )

base_volcano <- ggplot(brca1_kd, aes(x = logFC, y = -log10(P.Value))) + 
  labs(x = "log2(fold change)", y = "-log10(p-value)") + 
  theme_bw() + 
  theme(
    axis.ticks = element_line(colour = "black"),
    axis.text = element_text(colour = "black")
  )

# basic volcano plot
.plt <- base_volcano + 
  geom_point(aes(colour = change), size = 1) + 
  geom_hline(
    yintercept = -log10(0.05),
    linetype = "longdash",
    linewidth = 0.5,
    colour = "#7f7f7f"
  ) + 
  geom_vline(
    xintercept = c(-0.5, 0.5),
    linetype = "longdash",
    linewidth = 0.5,
    colour = "#7f7f7f"
  ) + 
  geom_text_repel(aes(label = sig_genes), max.overlaps = Inf) + 
  labs(subtitle = "BRCA1 knockdown volcano plot") + 
  scale_colour_manual(values = c(RColorBrewer::brewer.pal(3, "Set1"), "black")) + 
  theme(
    legend.position = "none"
  )
ggsave(
  "brca1_kd_volcano.png",
  plot = .plt,
  width = 4,
  height = 5,
  units = "in",
  dpi = 600
)

# gene label colour
gene_col <- RColorBrewer::brewer.pal(4, "Set1")[4]

# label differentiation markers
.plt <- base_volcano + 
  geom_point(
    data = brca1_kd %>% arrange(!is.na(diff_genes)),
    mapping = aes(colour = is.na(diff_genes)),
    size = 1
  ) + 
  geom_hline(
    yintercept = -log10(0.05),
    linetype = "longdash",
    linewidth = 0.5,
    colour = "#7f7f7f"
  ) + 
  geom_vline(
    xintercept = c(-0.5, 0.5),
    linetype = "longdash",
    linewidth = 0.5,
    colour = "#7f7f7f"
  ) + 
  geom_label_repel(aes(label = diff_genes), colour = gene_col, size = 5, max.overlaps = Inf) + 
  labs(subtitle = "Differentiation markers") + 
  scale_colour_manual(values = c(gene_col, "black")) +
  theme(legend.position = "none")
ggsave(
  "brca1_kd_diff.png",
  plot = .plt,
  width = 4,
  height = 5,
  units = "in",
  dpi = 600
)

# label regulon markers
.plt <- base_volcano + 
  geom_point(
    data = brca1_kd %>% arrange(!is.na(reg_genes)),
    mapping = aes(colour = is.na(reg_genes)),
    size = 1
  ) + 
  geom_hline(
    yintercept = -log10(0.05),
    linetype = "longdash",
    linewidth = 0.5,
    colour = "#7f7f7f"
  ) + 
  geom_vline(
    xintercept = c(-0.5, 0.5),
    linetype = "longdash",
    linewidth = 0.5,
    colour = "#7f7f7f"
  ) + 
  labs(subtitle = "Regulon target genes") + 
  scale_colour_manual(values = c(gene_col, "black")) + 
  theme(legend.position = "none")
ggsave(
  "brca1_kd_reg.png",
  plot = .plt,
  width = 4,
  height = 5,
  units = "in",
  dpi = 600
)
