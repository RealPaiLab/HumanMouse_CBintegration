# ==============================================================================
# These figures were generated for the OICR informatics lunch and learn on
# January 31, 2024. See also figures from `20240129-taylor_lab_meeting`.
# ==============================================================================

# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(tidyverse)
library(patchwork)
library(ggrepel)
library(ggpubr)

# load Cavalli data
cavalli_dir <- "/isilon/CBL_scRNAseq-archived/data/src/medulloblastoma-genomics/RNA/Cavalli_2017"
cavalli_meta <- readxl::read_excel(
  path = file.path(cavalli_dir, "mmc2.xlsx"),
  sheet = 1,
  skip = 1,
  col_names = TRUE
) %>% 
  select(1:18,48,49)
cavalli_expr <- read_delim(
  file = file.path(cavalli_dir, "GSE85217_M_exp_763_MB_SubtypeStudy_TaylorLab.txt.gz"),
  delim = "\t",
  col_names = TRUE
) %>% 
  filter(HGNC_symbol_from_ensemblv77 == "BRCA1") %>% 
  select(-c(1:3,5)) %>% 
  column_to_rownames(var = "HGNC_symbol_from_ensemblv77") %>% 
  t() %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "Study_ID")
  
# ------------------------------------------------------------------------------
# BRCA1 expression in i17q tumours (Group 4 ONLY as suggested by Michael
# Taylor); code from Xinghan

cavalli_meta <- mutate(
  cavalli_meta,
  i17q = case_when(
    `17p` == "deletion" & `17q` == "gain" ~ "yes",
    TRUE ~ "no"
  )
)
cavalli_meta$Subgroup <- factor(
  cavalli_meta$Subgroup,
  levels = c("WNT", "SHH", "Group3", "Group4")
)

cavalli_merge <- merge(cavalli_meta, cavalli_expr, by = "Study_ID")

.plt <- ggplot(filter(cavalli_merge, Subgroup %in% c("Group4")), aes(x = i17q, y = BRCA1)) +
# .plt <- ggplot(cavalli_merge, aes(x = i17q, y = BRCA1)) + 
  geom_violin() + 
  geom_jitter(width = 0.25, alpha=0.25) + 
  labs(y = "normalized log2(BRCA1 expression)") + 
  # facet_grid(cols = vars(Subgroup)) + 
  stat_compare_means(
    method = "wilcox.test",
    comparisons = list(c("yes", "no")),
    label = "p.format",
    tip.length = 0
  ) + 
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.08))) + 
  theme_classic() + 
  theme(
    axis.text = element_text(colour = "black"),
    axis.ticks = element_line(colour = "black")
  )

ggsave(
  filename = "brca1_expr_i17q.png",
  plot = .plt,
  width = 6,
  height = 4,
  units = "in",
  dpi = 600
)
