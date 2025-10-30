# ==============================================================================
# Differential gene expression fo MB tumour subtypes.
# ==============================================================================

library(argparse)
library(tidyverse)
library(Seurat)

# parse command line arguments
parser <- ArgumentParser()

parser$add_argument(
  # medulloblastoma Seurat object
  "--srat_rds",
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
    "--srat_rds", "/.mounts/labs/pailab/private/icheong/CBL_scRNAseq/results/tumour/Vladoiu/20230510/mb_mnn.rds",
    "--out_dir", "/.mounts/labs/pailab/private/icheong/CBL_scRNAseq/results/tumour/Vladoiu/20240917"
  ))
} else {
  arg_list <- parser$parse_args()
}

message(sprintf("***Saving files to %s***", arg_list$out_dir))
if (!dir.exists(arg_list$out_dir)) {
  dir.create(arg_list$out_dir, recursive = TRUE)
}

# ------------------------------------------------------------------------------
# load Seurat object

mb_srat <- readRDS(arg_list$srat_rds) %>%
  # recorrect SCT counts
  PrepSCTFindMarkers()

# set levels for MB subtypes
mb_srat$subtype <- factor(mb_srat$subtype, levels = c("SHH", "G3", "G4"))

# ------------------------------------------------------------------------------
# differential gene expression for MB clusters

# differential gene expression
message("***Running `FindAllMarkers` for MB clusters***")
diff_exp <- FindAllMarkers(
  mb_srat,
  assay = "SCT",
  logfc.threshold = 0,
  min.pct = 0.01
)

# save results
write_csv(
  x = diff_exp,
  file = file.path(arg_list$out_dir, "mb_cluster_diff_exp.csv")
)

# ------------------------------------------------------------------------------
# differential gene expression for MB subtypes

# set active idents to MB subtypes
Idents(mb_srat) <- "subtype"

# differential gene expression
message("***Running `FindAllMarkers` for MB subtypes***")
diff_exp <- FindAllMarkers(
  mb_srat,
  assay = "SCT",
  logfc.threshold = 0,
  min.pct = 0.01
)

# save results
write_csv(
  x = diff_exp,
  file = file.path(arg_list$out_dir, "mb_subtype_diff_exp.csv")
)

# ------------------------------------------------------------------------------

message("\n***SESSION INFO***\n")
print(sessionInfo())
