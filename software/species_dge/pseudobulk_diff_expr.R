# ==============================================================================
# Make pseudobulk and run differential gene expression across species.
# ==============================================================================

# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(argparse)
library(tidyverse)
library(ggrepel)
library(Seurat)
library(scuttle)
library(edgeR)

# parse command line arguments
parser <- ArgumentParser()

parser$add_argument(
  # human Seurat RDS file path
  "--human_srat",
  default = NULL,
  required = TRUE
)
parser$add_argument(
  # mouse Seurat RDS file path
  "--mouse_srat",
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
    "--human_srat", "/isilon/CBL_scRNAseq-archived/data/human/Aldinger/glutamatergic_dev_Liam.RDS",
    "--mouse_srat", "/CBL_scRNAseq/results/mouse/Vladoiu/merged_seurat_RLonly.rds",
    "--out_dir", "/CBL_scRNAseq/results/species_dge/20230829/"
  ))
} else {
  arg_list <- parser$parse_args()
}

message(sprintf("Saving files to %s", getwd()))
if (!dir.exists(arg_list$out_dir)) {
  dir.create(arg_list$out_dir, recursive = TRUE)
}

# load functions
source("/CBL_scRNAseq/software/mouse/Vladoiu/add_annotations.R")
source("/CBL_scRNAseq/software/utilities/convert_genes.R")
source("/CBL_scRNAseq/software/utilities/plotting.R")

# ------------------------------------------------------------------------------
# load data

human_srat <- read_rds(arg_list$human_srat)
mouse_srat <- read_rds(arg_list$mouse_srat) %>% 
  # custom function from `add_annotations.R`
  label_vladoiu_cells()

# ------------------------------------------------------------------------------
# convert mouse to orthologous human genes

orth_gene_list <- load_orth_genes("/CBL_scRNAseq/results/integrated/hgnc_mgi_orth_genes.csv")

mouse_srat <- rename_genes(mouse_srat, orth_gene_list$MGI.symbol, orth_gene_list$HGNC.symbol)

# ------------------------------------------------------------------------------
# merge human and mouse counts/metadata

comb_srat <- purrr::map2(
  .x = list(human_srat, mouse_srat),
  .y = c("human", "mouse"),
  .f = \(srat, species) {
    DefaultAssay(srat) <- "RNA"
    srat$species <- species
    srat <- DietSeurat(srat, counts = TRUE, assays = "RNA")
    return(srat)
  }
)
comb_srat <- merge(comb_srat[[1]], comb_srat[[2]]) %>% 
  # keep only orthologous genes
  subset(features = orth_gene_list$HGNC.symbol)

# convert back to factors
comb_srat$mouse_cell_type <- factor(
  comb_srat$mouse_cell_type,
  levels = levels(mouse_srat$mouse_cell_type)
)
comb_srat$new_cell_type <- factor(
  comb_srat$new_cell_type,
  levels = levels(human_srat$new_cell_type)
)
comb_srat$species <- factor(
  comb_srat$species,
  levels = c("human", "mouse")
)

# combine cell types into one columns
comb_srat@meta.data <- mutate(
  .data = comb_srat[[]],
  merge_cell_type = case_when(
    !is.na(new_cell_type) == TRUE ~ as.character(new_cell_type),
    TRUE ~ as.character(mouse_cell_type)
  )
)

# ------------------------------------------------------------------------------
# normalize + log transform

comb_srat <- NormalizeData(
  comb_srat,
  normalization.method = "LogNormalize",
  scale.factor = 10000
)


### DGE >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


#' Pseudobulk a Seurat object.
#'
#' @param srat Seurat object.
#' @param subset_expr Expression passed to `subset` the Seurat object. See
#'   Seurat documentation.
#' @param pb_group Metadata column(s) on which to pseudobulk the cells. For
#'   example, if you have a "species" column with human and mouse cells, then
#'   for each gene, it will add all the human counts together and all the mouse
#'   counts.
#'
#' @return A pseudobulked `SingleCellExperiment` object.
#' 
srat_to_pb <- function(
    srat,
    subset_expr,
    pb_group
) {
  # subset the Seurat object for the specific cell types
  subset_expr <- substitute(subset_expr)
  srat <- subset(srat, subset = !!subset_expr)
  
  # convert to SCE object
  sce <- as.SingleCellExperiment(srat)
  
  # pseudobulk
  groups <- colData(sce)[, pb_group, drop = FALSE]
  pb <- aggregateAcrossCells(sce, ids = groups)
  
  # add column names to pseudobulk
  if (length(pb_group) > 1) {
    df <- colData(pb)[, pb_group] %>% as.data.frame()
    cnames <- apply(df, MARGIN = 1, FUN = paste0, collapse = ".")
    colnames(pb) <- cnames
  } else {
    cnames <- colData(pb)[, pb_group]
    colnames(pb) <- cnames
  }
  
  return(pb)
}


#' Make `DGEList` object from pseudobulked `SingleCellExperiment` object.
#'
#' @param pb A `SingleCellExperiment` object with pseudobulk counts.
#'
#' @return A `DGEList` object.
#'
pb_to_dgelist <- function(
  pb  
) {
  dge_obj <- DGEList(
    counts = assay(pb),
    samples = colData(pb),
    group = rownames(colData(pb))
  )
  return(dge_obj)
}


#' Run `edgeR` differential gene expression using `exactTest`.
#'
#' @param dge_obj A `DGEList` object.
#' @param group Vector containing sample groups. Passed to `filterByExpr`.
#' @param bcv Biological coefficient of variation, default at 0.4. See `edgeR`
#'   documentation for more information. Briefly, a higher BCV means less
#'   significant *p*-values and FDRs.
#'
#' @return A `TopTags` object containing the FDRs of all genes tested.
#' 
run_dge <- function(
  dge_obj,
  group,
  bcv = 0.4
) {
  # filter out lowly expressed genes
  keep_genes <- filterByExpr(dge_obj, group = group)
  dge_obj <- dge_obj[keep_genes, ]
  message(sprintf(
    "Filtered out %s/%s genes with `filterByExpr`",
    sum(!keep_genes),
    length(keep_genes)
  ))
  
  # edgeR normalization of pseudobulk counts
  dge_obj <- calcNormFactors(dge_obj)
  
  message(sprintf("Testing %s vs. %s", group[2], group[1]))
  
  # differential expression
  extest <- exactTest(dge_obj, dispersion = bcv^2)
  results <- topTags(extest, n = Inf)
  
  return(results)
  
  
  # example of normal DGE (if you have replicates of samples)
  # rl_design <- model.matrix(~species, data = rl_dge$samples)
  # rl_dge <- estimateDisp(rl_dge, design = rl_design)
  # rl_ftest <- glmQLFit(rl_dge, design = rl_design, robust = TRUE) %>% 
  #   glmQLFTest()
  # rl_res <- topTags(rl_ftest, n = Inf)
}


#' Save `TopTags` object to `.csv` and `.rds` files.
#'
#' @param dge_res A `TopTags` object.
#' @param out_dir Output directory.
#' @param file_no_ext File name without the extension.
#'
#' @return None.
#'
save_dge_results <- function(
  dge_res,
  out_dir,
  file_no_ext
) {
  write_csv(
    x = dge_res %>% as.data.frame() %>% rownames_to_column(var = "gene"),
    file = file.path(out_dir, paste0(file_no_ext, ".csv"))
  )
  
  write_rds(
    x = dge_res,
    file = file.path(out_dir, paste0(file_no_ext, ".rds"))
  )
}


#' Add column named of gene labels for downstream volcano plot.
#'
#' @param dge_res A `TopTags` object.
#' @param gene_list List of genes you want to label.
#' @param top How many of the top genes should be labelled.
#'
#' @return A dataframe of the original object with a `gene_label` column.
#'
add_gene_labels <- function (
  dge_res,
  gene_list,
  top = 25
) {
  df <- as.data.frame(dge_res) %>% 
    rownames_to_column("gene") %>% 
    mutate(
      gene_label = case_when(
        gene %in% head(.$gene, top) ~ gene,
        gene %in% gene_list ~ gene,
        TRUE ~ ""
      )
    )
  return(df)
}


# gene list for volcano plot labelling
gene_list <- c("ATOH1",
               "BRCA1",
               "EOMES",
               "ERBB4",
               "LMX1A",
               "OTX2",
               "PAX6",
               "RBFOX3",
               "RELN")

# set BCV to 0.5 (dispersion = BCV^2)
bcv = 0.5


# ------------------------------------------------------------------------------
# DGE for RL-VZ

rlvz_pb <- srat_to_pb(
  srat = comb_srat,
  subset = new_cell_type == "RL-VZ" | 
    mouse_cell_type %in% c("Unipolar brush cell and GCP progenitor",
                           "Unipolar brush cell precursors"),
  pb_group = "species"
)

rlvz_dge <- pb_to_dgelist(rlvz_pb)

rlvz_res <- run_dge(rlvz_dge, group = colnames(rlvz_pb), bcv = bcv)

save_dge_results(rlvz_res, arg_list$out_dir, "01-rl_vz")


# make volcano plot
rlvz_df <- add_gene_labels(rlvz_res, gene_list = gene_list)

plt <- make_volcano(
  data = as.data.frame(rlvz_df),
  log_fc = logFC,
  log_pval = -log10(FDR),
  gene_labels = gene_label,
  log_pval_thresh = -log10(0.05)
) + 
  labs(title = "<-- human RL-VZ | mouse RL -->", y = "-log10 FDR")

ggsave(
  "01-rl_vz_volcano.png",
  plot = plt,
  path = arg_list$out_dir,
  width = 8,
  height = 8,
  units = "in",
  dpi = 600
)

# ------------------------------------------------------------------------------
# DGE for RL-SVZ

rlsvz_pb <- srat_to_pb(
  srat = comb_srat,
  subset = new_cell_type == "RL-SVZ" | 
    mouse_cell_type %in% c("Unipolar brush cell and GCP progenitor",
                           "Unipolar brush cell precursors"),
  pb_group = "species"
)

rlsvz_dge <- pb_to_dgelist(rlsvz_pb)

rlsvz_res <- run_dge(rlsvz_dge, group = colnames(rlsvz_pb), bcv = bcv)

save_dge_results(rlsvz_res, arg_list$out_dir, "02-rl_svz")


# make volcano plot
rlsvz_df <- add_gene_labels(rlsvz_res, gene_list = gene_list)

plt <- make_volcano(
  data = as.data.frame(rlsvz_df),
  log_fc = logFC,
  log_pval = -log10(FDR),
  gene_labels = gene_label,
  log_pval_thresh = -log10(0.05)
) + 
  labs(title = "<-- human RL-SVZ | mouse RL -->", y = "-log10 FDR")

ggsave(
  "02-rl_svz_volcano.png",
  plot = plt,
  path = arg_list$out_dir,
  width = 8,
  height = 8,
  units = "in",
  dpi = 600
)

# ------------------------------------------------------------------------------
# DGE for UBCs

ubc_pb <- srat_to_pb(
  srat = comb_srat,
  subset = new_cell_type %in% c("Early UBCs", "Late UBCs") | 
    mouse_cell_type == "Unipolar brush cells",
  pb_group = "species"
)

ubc_dge <- pb_to_dgelist(ubc_pb)

ubc_res <- run_dge(ubc_dge, group = colnames(ubc_pb), bcv = bcv)

save_dge_results(ubc_res, arg_list$out_dir, "03-ubc")


# make volcano plot
ubc_df <- add_gene_labels(ubc_res, gene_list = gene_list)

plt <- make_volcano(
  data = as.data.frame(ubc_df),
  log_fc = logFC,
  log_pval = -log10(FDR),
  gene_labels = gene_label,
  log_pval_thresh = -log10(0.05)
) + 
  labs(title = "<-- human UBC | mouse UBC -->", y = "-log10 FDR")

ggsave(
  "03-ubc_volcano.png",
  plot = plt,
  path = arg_list$out_dir,
  width = 8,
  height = 8,
  units = "in",
  dpi = 600
)

# ------------------------------------------------------------------------------
# testing DGE for UBCs vs GNs in humans only

testing_pb <- srat_to_pb(
  srat = comb_srat,
  subset = figure_clusters %in% c("04-GN", "05-eCN/UBC"),
  pb_group = "figure_clusters"
)

testing_dge <- pb_to_dgelist(testing_pb)

testing_res <- run_dge(testing_dge, group = colnames(testing_pb), bcv = 0.4)

save_dge_results(testing_res, arg_list$out_dir, "testing-gn_ubc")


# make volcano plot
testing_df <- add_gene_labels(testing_res, gene_list = gene_list)

plt <- make_volcano(
  data = as.data.frame(testing_df),
  log_fc = logFC,
  log_pval = -log10(PValue),
  gene_labels = gene_label,
  log_pval_thresh = -log10(0.05)
) + 
  labs(title = "<-- GN | UBC -->", y = "-log10 p-value")

ggsave(
  "testing-gn_ubc_volcano.png",
  plot = plt,
  path = arg_list$out_dir,
  width = 8,
  height = 8,
  units = "in",
  dpi = 600
)

# ------------------------------------------------------------------------------
# testing DGE for RL-SVZ vs UBCs in humans only

# combine Early and Late UBCs
comb_srat@meta.data <- mutate(
  comb_srat[[]],
  rl_ubcs = case_when(
    new_cell_type == "Early UBCs" ~ "UBCs",
    new_cell_type == "Late UBCs" ~ "UBCs",
    new_cell_type == "RL-SVZ" ~ "RL-SVZ",
    TRUE ~ NA_character_
  )
)

testing_pb <- srat_to_pb(
  srat = comb_srat,
  subset = (rl_ubcs %>% is.na %>% unname) == FALSE,# needs to be written this way, `!is.na(rl_ubcs)` gives an error
  pb_group = "rl_ubcs"
)

testing_dge <- pb_to_dgelist(testing_pb)

testing_res <- run_dge(testing_dge, group = colnames(testing_pb), bcv = 0.4)

save_dge_results(testing_res, arg_list$out_dir, "testing-rl_svz_ubc")


# make volcano plot
testing_df <- add_gene_labels(testing_res, gene_list = gene_list)

plt <- make_volcano(
  data = as.data.frame(testing_df),
  log_fc = logFC,
  log_pval = -log10(PValue),
  gene_labels = gene_label,
  log_pval_thresh = -log10(0.05)
) + 
  labs(title = "<-- RL-SVZ | UBC -->", y = "-log10 p-value")

ggsave(
  "testing-rl_svz_ubc_volcano.png",
  plot = plt,
  path = arg_list$out_dir,
  width = 8,
  height = 8,
  units = "in",
  dpi = 600
)

# ------------------------------------------------------------------------------

message("\n***SESSION INFO***\n")
print(sessionInfo())

