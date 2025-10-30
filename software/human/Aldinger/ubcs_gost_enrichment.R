# ==============================================================================
# Run GSEA with g:Profiler on genes upregulated/downregulated in human-specific
# UBCs to see which pathways are enriched.
# ==============================================================================

# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(argparse)
library(tidyverse)
library(ggrepel)
library(gprofiler2)

# parse command line arguments
parser <- ArgumentParser()

parser$add_argument(
  # differential expression file
  "--de_gene_file",
  default = NULL,
  required = TRUE
)
parser$add_argument(
  # output directory
  "--out_dir",
  default = NULL,
  required = TRUE
)

if (interactive()) {
  arg_list <- parser$parse_args(c(
    "--de_gene_file", "/CBL_scRNAseq/results/human/Aldinger/20230205/de_genes.csv",
    "--out_dir", "/CBL_scRNAseq/results/human/Aldinger/20240218"
  ))
} else {
  arg_list <- parser$parse_args()
}

message(sprintf("Saving files to %s", arg_list$out_dir))
if (!dir.exists(arg_list$out_dir)) {
  dir.create(arg_list$out_dir, recursive = TRUE)
}

# load functions
source("/CBL_scRNAseq/software/utilities/gprofiler2_helpers.R")
source("/CBL_scRNAseq/software/utilities/plotting.R")

# ------------------------------------------------------------------------------
# load differentially expressed genes in common vs. human-specific UBCs

de_genes <- read_csv(arg_list$de_gene_file)

human_ubc <- de_genes %>% 
  filter(p_val_adj < 0.05 & avg_log2FC < 0) %>% 
  pull(gene)
write_lines(
  human_ubc,
  file = file.path(arg_list$out_dir, "human_specific_ubc_genes.txt")
)

common_ubc <- de_genes %>% 
  filter(p_val_adj < 0.05 & avg_log2FC > 0) %>% 
  pull(gene)
write_lines(
  common_ubc,
  file = file.path(arg_list$out_dir, "common_ubc_genes.txt")
)

# re-make volcano plot
.plt <- de_genes %>% 
  mutate(
    gene_label = case_when(
      gene %in% head(.$gene, 30) ~ gene,
      TRUE ~ ""
    )
  ) %>% 
  make_volcano(
    data = .,
    log_fc = avg_log2FC,
    log_pval = -log10(p_val_adj),
    direction = NULL,
    gene_labels = gene_label,
    log_pval_thresh = -log10(0.05)
  ) + 
  labs(
    title = "<-- human-specific UBCs | common UBCs -->",
    x = "log2(fold change)",
    y = "-log10(adj. p-value)"
  )
ggsave(
  filename = "volcano.png",
  plot = .plt,
  path = arg_list$out_dir,
  width = 6,
  height = 6,
  units = "in",
  dpi = 600
)

# ------------------------------------------------------------------------------
# find enriched pathways

pwalk(
  .l = list(
    list(human_ubc, common_ubc),
    list("human_specific_ubc", "common_ubc")
  ),
  .f = \(gene_list, filename) {
    message(sprintf("***Running pathway enrichment analysis for %s***", filename))
    gost_res <- run_gost(
      query = gene_list,
      custom_bg = de_genes$gene,
      filename = file.path(arg_list$out_dir, filename)
    )
    
    gost_res2gem(gost_res = gost_res, phenotype = "+1") %>% 
      write_gem(file = paste0(
        file.path(arg_list$out_dir, filename),
        ".gem.txt"
      ))
  }
)

# ------------------------------------------------------------------------------

message("\n***SESSION INFO***\n")
print(sessionInfo())

