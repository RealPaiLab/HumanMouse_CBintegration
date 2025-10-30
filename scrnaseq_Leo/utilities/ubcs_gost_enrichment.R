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
  # background file
  "--background_file",
  default = NULL,
  required = FALSE
)
parser$add_argument(
  # output directory
  "--out_dir",
  default = NULL,
  required = TRUE
)

if (interactive()) {
  arg_list <- parser$parse_args(c(
    "--de_gene_file", "/.mounts/labs/pailab/private/llau/results/integrated/20240715/sig_ubc_markers.csv",
    "--background_file", "/.mounts/labs/pailab/private/llau/results/integrated/20240715/background.txt",
    "--out_dir", "/.mounts/labs/pailab/private/llau/results/integrated/20240711"
  ))
} else {
  arg_list <- parser$parse_args()
}

if (!dir.exists(arg_list$out_dir)) {
  dir.create(arg_list$out_dir, recursive = TRUE)
}

log_file <- file(paste0(arg_list$out_dir, "/log.txt"), open = "wt")

sink(log_file)
sink(log_file, type = "message")

message(sprintf("Saving files to %s", arg_list$out_dir))
# load functions
source("/u/llau/software/mb_scrnaseq/MB_scRNAseq/software/utilities/gprofiler2_helpers.R")

# ------------------------------------------------------------------------------
# load differentially expressed genes in specified group

de_genes_list <- read_csv(arg_list$de_gene_file) %>%
  split(.$cluster)

if(!is.null(arg_list$background_file))
{
  message(paste0("Using provided background: ", arg_list$background_file))
  custom_bg <- read.csv(arg_list$background_file) %>%
    split(.$cluster)
} 

for(cluster in names(de_genes_list)){
  out_path <- paste0(arg_list$out_dir, "/", cluster)
  if (!dir.exists(out_path)) {
    dir.create(out_path, recursive = TRUE)
  }

  de_genes <- de_genes_list[[cluster]]

  if(!is.null(arg_list$background_file))
  {
    message(paste0("In cluster ", cluster, " there were ", nrow(custom_bg[[cluster]]), " tested genes"))
  } else {
    message(paste0("In cluster ", cluster, " there were ", nrow(de_genes), " tested genes"))
  }

  upregulated <- de_genes %>% 
    filter(p_val_adj < 0.05 & avg_log2FC > 0) %>% 
    pull(feature)
  write_lines(
    upregulated,
    file = file.path(out_path, paste0(cluster, "_upregulated_genes.txt"))
  )
  message(paste0("There are ", length(upregulated), " significant upregulated genes"))

  downregulated <- de_genes %>% 
    filter(p_val_adj < 0.05 & avg_log2FC < 0) %>% 
    pull(feature)
  write_lines(
    downregulated,
    file = file.path(out_path, paste0(cluster, "_downregulated_genes.txt"))
  )
  message(paste0("There are ", length(downregulated), " significant downregulated genes"))

  # ------------------------------------------------------------------------------
  # find enriched pathways

  pwalk(
    .l = list(
      list(upregulated, downregulated),
      list("upregulated", "downregulated")
    ),
    .f = \(gene_list, filename) {
      message(sprintf("***Running pathway enrichment analysis for %s***", filename))

      if(!is.null(arg_list$background_file))
      {
        bg <- custom_bg[[cluster]]$feature
      } else {
        bg <- de_genes$feature
      }

      gost_res <- run_gost(
        query = gene_list,
        significant = FALSE,
        organism = "gp__OUFl_gpIE_RyE",
        custom_bg = bg,
        filename = file.path(out_path, filename)
      )

      gost_res2gem(gost_res = gost_res, phenotype = "+1") %>% 
        write_gem(file = paste0(
          file.path(out_path, filename),
          ".gem.txt"
        ))
    }
  )
}
# ------------------------------------------------------------------------------

message("\n***SESSION INFO***\n")
print(sessionInfo())

sink(type = "message")
sink()