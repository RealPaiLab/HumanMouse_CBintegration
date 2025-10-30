library(argparse)
library(Seurat)
library(qs)
library(SingleR)
library(SingleCellExperiment)
library(scran)

# parse command line arguments
parser <- ArgumentParser()

parser$add_argument(
  # path to reference SCE object
  "--ref_sce",
  default = NULL,
  required = TRUE
)
parser$add_argument(
  # metadata column name for reference dataset cell types
  "--col_name",
  default = NULL,
  required = FALSE
)
parser$add_argument(
  # path to tumour SCE object
  "--tumour_sce",
  default = NULL,
  required = FALSE
)
parser$add_argument(
  # output directory
  "--out_dir",
  default = NULL,
  required = TRUE
)

arg_list <- parser$parse_args()

if (!dir.exists(arg_list$out_dir)) {
  dir.create(arg_list$out_dir, recursive = TRUE)
}

rl_sce <- qread(arg_list$ref_sce)
tumour_sce <- qread(arg_list$tumour_sce)

trained_model <- trainSingleR(
  rl_sce, 
  labels=rl_sce[[arg_list$col_name]], 
  de.method = "wilcox", 
  assay.type = "logcounts",
  BPPARAM=MulticoreParam(8)
)

predict <- classifySingleR(
  test = tumour_sce,
  trained = trained_model,
  assay.type = "logcounts",
  BPPARAM=MulticoreParam(8)
)


qsave(rl_sce, file.path(out_directory, "rl_sce.qs"))
qsave(tumour_sce, file.path(out_directory, "tumour_sce.qs"))
qsave(predict, file.path(out_directory, "predictions.qs"))