# ==============================================================================
# Calculate pseudoages of each cell and Manhattan distance between ages to
# determine the correspondence of developmental stages across mouse and human
# (see methods from Sepp et al. 2024).
# ==============================================================================

# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(argparse)
library(tidyverse)
library(Seurat)

# parse command line arguments
parser <- ArgumentParser()

parser$add_argument(
  # integrated Seurat object
  "--srat_qs",
  default = NULL,
  required = TRUE
)
parser$add_argument(
  # dataset name to use for reference
  "--ref_dataset",
  default = "Sepp_RL_mouse"
)
parser$add_argument(
  # output directory
  "--out_dir",
  default = NULL,
  required = TRUE
)
parser$add_argument(
  # calculate manhattan distance
  "--save_dist_mat",
  action = "store_true"
)

if (interactive()) {
  arg_list <- parser$parse_args(c(
    "--srat_qs", "/.mounts/labs/pailab/private/llau/results/integrated/20240516/fc_subset_analysis/20240514_cca_integ_subset.qs",
    "--ref_dataset", "Sepp_full_cerebellum_mouse",
    "--out_dir", "/.mounts/labs/pailab/private/icheong/CBL_scRNAseq/results/integrated/20240521",
    "--save_dist_mat"
  ))
} else {
  arg_list <- parser$parse_args()
}

message(sprintf("***Saving files to %s***", arg_list$out_dir))
if (!dir.exists(arg_list$out_dir)) {
  dir.create(arg_list$out_dir, recursive = TRUE)
}

# ------------------------------------------------------------------------------
# functions

source("/.mounts/labs/pailab/private/icheong/CBL_scRNAseq/software/utilities/pseudoage_helpers.R")

# ------------------------------------------------------------------------------
# load integrated Seurat object

srat <- qs::qread(arg_list$srat_qs)

# collect sample timepoints into two columns (one human, one mouse)
# `dataset_name` may be in the format of "<author>_RL_<species>" or
# "<author>_full_cerebellum_<species>", so use regex to check both options
srat@meta.data <- srat@meta.data %>% 
  mutate(
    # timepoints from human datasets
    human_age = case_when(
      str_detect(dataset_name, "^Aldinger\\w*human$") ~ age,
      str_detect(dataset_name, "^Luo\\w*human$") ~ paste(
        str_remove(string = batch, pattern = "PCW"),
        "PCW"
      ),
      str_detect(dataset_name, "^Sepp\\w*human$") ~ paste(
        str_remove(string = Stage, pattern = " wpc"),
        "PCW"
      )
    ),
    # timepoints from mouse datasets
    mouse_age = case_when(
      str_detect(dataset_name, "^Vladoiu\\w*mouse$") ~ str_remove(
        string = orig.ident,
        pattern = "Vladoiu-"
      ),
      str_detect(dataset_name, "^Sepp\\w*mouse$") ~ Stage
    )
  )

# convert ages to factor
srat$human_age <- fct(
  x = srat$human_age,
  levels = str_sort(unique(srat$human_age), numeric = TRUE, na_last = NA)
)
srat$mouse_age <- fct(
  x = srat$mouse_age,
  levels = str_sort(unique(srat$mouse_age), numeric = TRUE, na_last = NA)
)

# ------------------------------------------------------------------------------
# sample dataset to generate reference (defaults to Sepp mouse dataset)

# subset cells
sepp_ref <- subset(srat, dataset_name == arg_list$ref_dataset)

# drop empty levels in mouse_age column
sepp_ref$mouse_age <- fct_drop(sepp_ref$mouse_age)

# set index of stages
sepp_ref$mouse_age_idx <- as.numeric(sepp_ref$mouse_age)

# sample cell IDs by index and cell type
n_max <- 250 # maximum number of cells per cell type to sample
sepp_ref$index_cell_type <- paste0(
  sepp_ref$mouse_age_idx,
  "//",
  sepp_ref$common_cell_name
)
sepp_ref_cells <- purrr::map(
  .x = unique(sepp_ref$index_cell_type),
  .f = function(x) {
    cells <- rownames(sepp_ref[[]][sepp_ref$index_cell_type == x, ])
    set.seed(42)
    sample(x = cells, size = min(length(cells), n_max))
  }
) %>%
  unlist()
message(sprintf(
  "***Sampled %s/%s Sepp mouse cells to generate the reference timepoints***",
  length(sepp_ref_cells),
  nrow(sepp_ref[[]])
))

# select reference cells by cell IDs
sepp_ref <- subset(sepp_ref, cells = sepp_ref_cells)

# save reference cell metadata
message("***Saving Sepp mouse reference cell type developmental index to `sepp_mouse_ref.csv`***")
sepp_ref[[]] %>%
  rownames_to_column("cell_id") %>%
  select(cell_id, mouse_age, mouse_age_idx, index_cell_type) %>%
  write_csv(file = file.path(arg_list$out_dir, "sepp_mouse_ref.csv"))

# ------------------------------------------------------------------------------
# calculate pseudoage for each cell

# smaller samples for testing
if (FALSE) {
  sub_srat <- subset(srat, cells = sample(rownames(srat[[]]), 2500))
  sub_ref <- subset(sepp_ref, cells = sample(sepp_ref_cells, 250))
  pseudoage <- get_pseudostage(
    srat = sub_srat,
    ref_srat = sub_ref,
    ref_stages = "mouse_age_idx",
    reduction = "pca",
    k = 25,
    metric = "cosine"
  )
}

# set number of nearest neighbours to search for
k <- 25

# pseudoage_helpers.R
pseudoage <- get_pseudostage(
  srat = srat,
  ref_srat = sepp_ref,
  ref_stages = "mouse_age_idx",
  reduction = "pca",
  k = k,
  metric = "cosine"
)

# get actual human and mouse age
true_age <- srat[[]] %>%
  select(mouse_age, human_age)

# merge actual ages with pseudoages
all_age <- merge(
  x = true_age,
  y = as.data.frame(pseudoage),
  by = "row.names",
  all = TRUE
) %>%
  rename(cell_id = "Row.names")

# ------------------------------------------------------------------------------
# bin pseudoage and calculate Manhattan distance

message("***Saving pseudoages to `pseudoage.csv`***")
write_csv(
  x = all_age,
  file = file.path(arg_list$out_dir, "pseudoage.csv")
)

# calculate Manhattan distance for all pairwise developmental stages
if (arg_list$save_dist_mat) {
  md <- merge(
    x = srat[[]],
    y = select(.data = all_age, cell_id, pseudoage) %>% column_to_rownames("cell_id"),
    by = "row.names",
    all = TRUE
  ) %>%
    column_to_rownames("Row.names")

  # proportion of cells in each pseudoage bin per developmental timepoint
  message("***Binning pseudoages***")
  prop_bin <- bin_pseudoage(md, bin_col = "pseudoage")

  # Manhattan distance
  pseudoage_dist <- prop_bin %>%
    dist(method = "manhattan", diag = TRUE, upper = TRUE) %>%
    as.matrix()

  message("***Saving Manhattan distance matrix to `pseudoage_dist.csv`***")
  write.csv(
    pseudoage_dist,
    file = file.path(arg_list$out_dir, "pseudoage_dist.csv")
  )
}

# ------------------------------------------------------------------------------

message("\n***SESSION INFO***\n")
print(sessionInfo())
