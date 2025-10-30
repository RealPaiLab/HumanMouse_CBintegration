library(argparse)
library(tidyverse)
library(Seurat)
library(qs)

# Parse command line arguments
parser <- ArgumentParser()

parser$add_argument(
  # Input file
  "-file",
  action = "store",
  required = TRUE,
  help = 'Path to integration file to be subsetted'
)

parser$add_argument(
  # Subset size
  "-size",
  action = "store",
  required = TRUE,
  help = 'A non-negative integer giving the number of items to choose'
)

parser$add_argument(
  # Output directory
  "-output",
  action = "store",
  required = TRUE,
  help = 'Path to output directory'
)

arg_list <- parser$parse_args()

integ_file_path <- as.character(arg_list$file)
sample_size <- as.integer(arg_list$size)
out_directory <- as.character(arg_list$output)

# Setting the seed
set.seed(1) 

# Loading in the integration file
srat_integ <- qread(integ_file_path)

# Sampling and subsetting of integration file
sampled_cells <- sample(rownames(srat_integ@meta.data), sample_size, replace = FALSE, prob = NULL)
sampled_srat <- subset(srat_integ, cells = sampled_cells)

# Saving subsetted file
filename <- paste0(out_directory, "/", format(Sys.Date(), "%Y%m%d"), "_", sample_size, "_subset.qs")
qsave(sampled_srat, filename)