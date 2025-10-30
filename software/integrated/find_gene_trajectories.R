# ==============================================================================
# Determine gene expression trajectories to identify conserved and non-conserved
# genes (see methods from Sepp et al. 2024).
# ==============================================================================

library(argparse)
library(tidyverse)
library(Seurat)
library(Mfuzz)

# parse command line arguments
parser <- ArgumentParser()

parser$add_argument(
  # integrated Seurat object
  "--srat_qs",
  default = NULL,
  required = TRUE
)
parser$add_argument(
  # pseudotime values
  "--pseudotime_csv",
  default = NULL,
  required = TRUE
)
parser$add_argument(
  # which lineage to use (UBC or GCs?)
  "--lineage",
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
    "--srat_qs", "/.mounts/labs/pailab/private/llau/results/integrated/20240524/25_pc_without_luo/25_pc_rl.qs",
    "--pseudotime_csv", "/.mounts/labs/pailab/private/icheong/CBL_scRNAseq/results/integrated/20240611/integ_pseudotime/pseudotime_values.csv",
    "--lineage", "ubc",
    "--out_dir", "/.mounts/labs/pailab/private/icheong/CBL_scRNAseq/results/integrated/20240711"
  ))
} else {
  arg_list <- parser$parse_args()
}

message(sprintf("***Saving files to %s***", arg_list$out_dir))
if (!dir.exists(arg_list$out_dir)) {
  dir.create(arg_list$out_dir, recursive = TRUE)
}

# set number of pseudotime bins
n_bins <- 10

# set minimum number of cells in each pseudobulk group
min_cells_per_pseudobulk <- 25

# set cutoff for the Poisson variance mean ratio
pvmr_cutoff <- 1

# set number of fuzzy clusters to create
mfuzz_nclusters <- 4

# ------------------------------------------------------------------------------
# functions

source("/.mounts/labs/pailab/private/icheong/CBL_scRNAseq/software/utilities/cell_labelling.R")

# ------------------------------------------------------------------------------

# steps
# for each lineage of cells:
# - bin pseudotime
# - split by species
# - pseudobulk per replicate per bin (use dataset for replicates?)
# - merge pseudobulk counts in each bin (take the mean)
# - normalize to counts per million for each species with orthologous genes
# - find overdispersed genes (aka dynamic genes) in each species
# - keep only intersection of overdispersed genes across species
# - combine CPM across species into a matrix with genes in rows and pseudotime bins in columns (genes appear once for each species)
# - clustering with Mfuzz and determining (adjusted) p-value to determine similarity of trajectory of two genes
# - determine which genes are preserved vs. diverged

# ------------------------------------------------------------------------------
# load Seurat objects

# load integrated Seurat object
message(sprintf("***Loading in Seurat object at %s***", arg_list$srat_qs))
srat <- qs::qread(arg_list$srat_qs)

# save the orthologous features (genes) used for integration
orth_genes <- VariableFeatures(srat, assay = "integrated")

# set the resolution as the seurat clusters
srat$seurat_clusters <- srat$snn_res.0.4

# clean up the dataset column
srat$dataset_name <- str_remove(
  string = srat$dataset_name,
  pattern = "full_cerebellum_"
)

# add annotations for cell type and UBC subclusters (see `cell_labelling.R`)
srat <- label_rl_lineage_integration(srat)

# subset by lineage (UBC or GC)
if (arg_list$lineage == "ubc") {
  # get cells in the UBC lineage
  message("***Determining gene expression trajectories for the UBC lineage***")
  cell_types <- c("RL", paste0("UBC_", 0:5))
} else if (arg_list$lineage == "gc") {
  # get cells in the granule cell lineage
  message("***Determining gene expression trajectories for the granule cell lineage***")
  cell_types <- c("RL", "GCP", "GC")
} else {
  stop("`--lineage` must be one of `ubc` or `gc`")
}
cell_ids <- rownames(srat[[]][srat$cell_type_annot %in% cell_types, ])
srat <- srat[, cell_ids]

# ------------------------------------------------------------------------------
# load pseudotime and make bins

# load pseudotime values
pseudotime <- read_csv(arg_list$pseudotime_csv) %>%
  column_to_rownames("cell_id")

# the control cells (oligodendrocytes, microglia, endothelial cells) were
# excluded from the pseudotime analysis; keep only the cells that were included
# in the pseudotime analysis
srat <- srat[, rownames(pseudotime)]

# add pseudotime values to the Seurat metadata
md <- merge(
  x = srat[[]],
  y = pseudotime,
  by = "row.names",
  all.x = TRUE
) %>%
  column_to_rownames("Row.names")

# now bin the pseudotime values
message(sprintf("***Binning pseudotime values into %s***", n_bins))
md$pseudotime_bins <- cut(x = md$pseudotime, breaks = n_bins, labels = FALSE)
# md$pseudotime_bins <- cut_number(x = md$pseudotime, n = n_bins, labels = FALSE)

# sort new metadata to original order and add it back to Seurat object
md <- md[rownames(srat[[]]), ]
srat@meta.data <- md

# ------------------------------------------------------------------------------
# split by species and pseudobulk

# keep only orthologous genes that were used in the integration
srat <- srat[orth_genes, ]

# separate human and mouse
srat_list <- SplitObject(srat, split.by = "species")

# create pseudobulk per pseudotime bin (using datasets as replicates)
pb_reps <- map(
  .x = srat_list,
  .f = \(s) {
    # create group for each dataset and pseudotime bin
    s$grp <- paste(s$dataset_name, s$pseudotime_bins, sep = "//")

    # remove pseudobulks with few cells
    grp_table <- table(s$grp)
    keep_grp <- names(grp_table)[grp_table >= min_cells_per_pseudobulk]
    s <- s[, rownames(s[[]][s$grp %in% keep_grp, ])]

    # calculate pseudobulk for each dataset/bin
    pb <- AggregateExpression(
      s,
      assays = "RNA",
      group.by = "grp",
      slot = "counts"
    ) %>%
      pluck("RNA")

    colnames(pb) <- unname(colnames(pb))

    return(pb)
  }
)

# merge pseudobulks in same pseudotime bin (take the average)
pb_avg <- map(
  .x = pb_reps,
  .f = \(pb) {
    # group by pseudotime bin
    grp <- colnames(pb) %>%
      str_split_fixed(string = ., pattern = "//", n = 2) %>%
      .[, 2]

    # calculate mean expression in each bin
    pb <- tapply(
      X = seq_len(ncol(pb)),
      INDEX = grp,
      FUN = \(i) rowMeans(pb[, i, drop = FALSE]),
      simplify = FALSE
    ) %>%
      do.call(cbind, .)

    return(pb)
  }
)

# keep only pseudotime bins that are present in both species
keep_bins <- map(pb_avg, colnames) %>% purrr::reduce(intersect)
pb_avg <- map(
  .x = pb_avg,
  .f = \(pb) pb[, keep_bins]
)
message(sprintf(
  "***Using %s of %s pseudotime bins for downstream analysis***",
  length(keep_bins),
  n_bins
))

# normalize pseudobulked genes to counts per million (CPM)
pb_norm <- map(
  .x = pb_avg,
  .f = \(pb) {
    # normalize by counts in each bin and multiply by 1e6
    t(t(pb) / colSums(pb)) * 1e6
  }
)

# ------------------------------------------------------------------------------
# find highly variable (overdispersed) genes

# calculate variance of each gene
gene_vars <- map(
  .x = pb_avg,
  .f = \(pb) {
    # normalize by counts in each bin
    t(t(pb) / colSums(pb)) %>%
      # get variance of normalized counts for each gene
      matrixStats::rowVars()
  }
)

# calculate mean of each gene
gene_means <- map(
  .x = pb_avg,
  .f = \(pb) {
    # normalize by counts in each bin
    t(t(pb) / colSums(pb)) %>%
      # get mean of normalized counts for each gene
      base::rowMeans()
  }
)

# calculate expected variance to mean ratio (pVMR)
gene_pvmr <- map(
  .x = pb_avg,
  .f = \(pb) mean(1 / colSums(pb))
)

# calculate center of mass of each gene (on average, which pseudotime bin each
# gene is "most expressed" in)
gene_com <- map(
  .x = pb_norm,
  .f = \(pb) {
    # the code from the paper uses Matrix::rowSums instead of base::rowSums - no
    # idea why though...
    colSums(t(pb) * seq_len(ncol(pb))) / Matrix::rowSums(pb)
  }
)

# calculate pseudotime bin specificity of gene expression
gene_time_tau <- map(
  .x = pb_norm,
  .f = \(pb) {
    rowSums(1 - (pb / matrixStats::rowMaxs(pb))) / (ncol(pb) - 1)
  }
)

# calculate number of cells expressing each gene (raw counts > 0)
gene_num_cells <- map(
  .x = srat_list,
  .f = \(s) rowSums(s[["RNA"]]@counts > 0)
)

# table of gene stats
gene_stats <- pmap(
  .l = list(
    gene_means,
    gene_vars,
    gene_pvmr,
    gene_num_cells,
    names(gene_num_cells),
    gene_com,
    gene_time_tau
  ),
  .f = \(
    gene_means,
    gene_vars,
    gene_pvmr,
    gene_num_cells,
    species,
    gene_com,
    gene_time_tau
  ) {
    tibble(
      gene_name = names(gene_means),
      mean = gene_means,
      var = gene_vars,
      centre_of_mass = gene_com,
      pvmr = gene_pvmr,
      time_tau = gene_time_tau,
      num_cells = gene_num_cells,
      species = species
    )
  }
) %>%
  do.call(what = rbind, args = .) %>%
  # calculate variance mean ratio for each gene
  mutate(vmr = var / mean) %>%
  # determine if the VMR is greater than expected
  mutate(select_gene = vmr > pvmr * pvmr_cutoff)

# filter for highly variable (overdispersed) genes
keep_genes <- gene_stats %>%
  dplyr::filter(select_gene) %>%
  group_by(gene_name) %>%
  filter(n() == 2) %>%
  pull(gene_name) %>%
  unique()
message(sprintf(
  "***Keeping %s overdispersed genes out of %s total genes***",
  length(keep_genes),
  length(orth_genes)
))

# scale the CPM for the highly variable genes
pb_scaled <- map(
  .x = pb_norm,
  .f = \(pb) t(scale(t(pb[keep_genes, ])))
)

# save gene statistics and highly variable genes
write_csv(
  x = gene_stats,
  file = file.path(arg_list$out_dir, "gene_stats.csv")
)

# save scaled pseudobulked gene expression (one list element for each species)
saveRDS(
  object = pb_scaled,
  file = file.path(arg_list$out_dir, "scaled_pseudobulk_mat_list.rds")
)

# ------------------------------------------------------------------------------
# soft clustering with Mfuzz

# create matrix/object for Mfuzz clustering
eset <- map2(
  .x = pb_scaled,
  .y = names(pb_scaled),
  .f = \(pb, species) {
    rownames(pb) <- paste(species, rownames(pb), sep = "_")
    return(pb)
  }
) %>%
  do.call(what = rbind, args = .) %>%
  ExpressionSet()

# estimate optimal fuzzifier parameter
estimated_m <- mestimate(eset)
message(sprintf("***Estimated optimal fuzzifier m = %.4f***", estimated_m))

# run Mfuzz
set.seed(42)
mf <- mfuzz(
  eset = eset,
  centers = mfuzz_nclusters,
  m = estimated_m
)

# plot expression vs. pseudotime for each of the clusters
num_col <- ceiling(sqrt(mfuzz_nclusters))
num_row <- ceiling(mfuzz_nclusters / num_col)
png(
  filename = file.path(arg_list$out_dir, "mfuzz_plot.png"),
  width = 4 * num_col,
  height = 4 * num_row,
  units = "in",
  res = 300
)
mfuzz.plot(
  eset = eset,
  cl = mf,
  mfrow = c(num_row, num_col),
  new.window = FALSE
)
par(mfrow = c(1, 1)) # reset plotting
dev.off()

# isolate cluster values (scores)
mf_mem_scores <- as_tibble(mf$membership)
colnames(mf_mem_scores) <- paste0("mf_mem_", colnames(mf_mem_scores))

# get optimal cluster assignments
mf_df <- tibble(
  species_gene_name = names(mf$cluster),
  mf_cluster = mf$cluster
) %>%
  bind_cols(mf_mem_scores) %>%
  separate_wider_delim(
    cols = species_gene_name,
    delim = "_",
    names = c("species", "gene")
  )

# calculate p-value for similarity of trajectory between species
mf_pval <- mf_df %>%
  pivot_longer(
    cols = starts_with("mf_mem_"),
    names_to = "mf_mem",
    values_to = "mf_mem_score"
  ) %>%
  # filter for genes with max score > 0.5 in both species
  group_by(gene, species) %>%
  filter(max(mf_mem_score) > 0.5) %>%
  group_by(gene) %>%
  filter(length(unique(species)) == 2) %>%
  # remove column
  select(-mf_cluster) %>%
  # calculate p-value and FDR for each gene
  pivot_wider(names_from = species, values_from = mf_mem_score) %>%
  summarise(
    pval = sum(human * mouse),
    mf_clust_h = max(human),
    mf_clust_m = max(mouse)
  ) %>%
  mutate(
    fdr = p.adjust(p = pval, method = "BH")
  )

# classify orthologs
diff_fdr <- 0.05 # FDR cutoff
same_fdr <- 0.5 # similarity cutoff
mf_pval <- mf_pval %>%
  left_join(
    y = mf_df %>%
      select(gene, species, mf_cluster) %>%
      pivot_wider(names_from = species, values_from = mf_cluster) %>%
      dplyr::rename(cluster_human = human, cluster_mouse = mouse),
    by = join_by(gene)
  ) %>%
  mutate(mf_class = case_when(
    fdr < diff_fdr ~ "diverged",
    fdr > same_fdr & cluster_human == cluster_mouse ~ "conserved",
    .default = "not assigned"
  ))

# save ExpressionSet object
saveRDS(
  object = eset,
  file = file.path(arg_list$out_dir, "scaled_pseudobulk_eset.rds")
)

# save Mfuzz clustering output
saveRDS(
  object = mf,
  file = file.path(arg_list$out_dir, "mfuzz_raw_output.rds")
)

# save Mfuzz clustering score table
write_csv(
  x = mf_df,
  file = file.path(arg_list$out_dir, "mfuzz_membership_scores.csv")
)

# save Mfuzz clustering p-values and assignment
write_csv(
  x = mf_pval,
  file = file.path(arg_list$out_dir, "mfuzz_clustering_pval.csv")
)

# ------------------------------------------------------------------------------

message("\n***SESSION INFO***\n")
print(sessionInfo())
