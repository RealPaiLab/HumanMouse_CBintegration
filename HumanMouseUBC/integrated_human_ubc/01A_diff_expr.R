# ==============================================================================
# Run differential gene expression on integrated human UBC subclusters.
# ==============================================================================

library(argparse)
library(tidyverse)
library(patchwork)
library(Seurat)
library(ComplexHeatmap)

# parse command line arguments
parser <- ArgumentParser()

parser$add_argument(
  # Seurat file
  "--srat_file",
  default = NULL,
  required = TRUE
)
parser$add_argument(
  # differential expression test (see Seurat FindMarkers for list of tests)
  "--test_use",
  default = "wilcox",
  required = FALSE
)
parser$add_argument(
  # latent vars to include in model (see Seurat FindMarkers)
  "--latent_vars",
  action = "extend",
  nargs = "*",
  default = NULL,
  required = FALSE
)
parser$add_argument(
  # if true, run `FindConservedMarkers`
  "--conserved",
  action = "store_true"
)
parser$add_argument(
  # output directory
  "--out_dir",
  default = NULL,
  required = TRUE
)

if (interactive()) {
  # for testing and troubleshooting
  arg_list <- parser$parse_args(c(
    "--srat_file", "/.mounts/labs/pailab/public/HumanMouseUBC/data/UBC.Harmony.RDS",
    "--test_use", "MAST",
    "--latent_vars", "sex", "dataset_name",
    "--out_dir", "/.mounts/labs/pailab/private/projects/HumanMouseUBC/integrated_human_ubc/20250331/"
  ))
} else {
  arg_list <- parser$parse_args()
}

message(sprintf("Saving files to %s", arg_list$out_dir))
if (!dir.exists(arg_list$out_dir)) {
  dir.create(arg_list$out_dir, recursive = TRUE)
}

# ------------------------------------------------------------------------------
# functions

source("/.mounts/labs/pailab/private/icheong/CBL_scRNAseq/software/utilities/plotting.R")

run_all_de <- function(
  srat,
  clusters,
  conserved = arg_list$conserved,
  test_use = arg_list$test_use,
  latent_vars = arg_list$latent_vars
) {
  if (conserved) {
    message(sprintf(
      "***Using `FindConservedMarkers` to run differential gene expression with `test.use = %s` and `latent.vars = %s`***",
      test_use,
      latent_vars
    ))
    de_func <- FindConservedMarkers
  } else {
    message(sprintf(
      "***Using `FindMarkers` to run differential gene expression with `test.use = %s` and `latent.vars = c(%s)`***",
      test_use,
      paste0(latent_vars, collapse = ", ")
    ))
    de_func <- FindMarkers
  }

  # run differential gene expression (each cluster vs. all other clusters)
  de_genes <- map(
    .x = clusters,
    .f = \(clust) {
      message(sprintf("***%s, n = %s cells***", clust, sum(Idents(srat) == clust)))

      # differential expression
      de_genes <- de_func(
        object = srat,
        ident.1 = clust,
        grouping.var = "dataset_name",
        assay = "SCT",
        logfc.threshold = 0,
        min.pct = 0.01,
        test.use = test_use,
        latent.vars = latent_vars
      ) %>%
        # add columns for cluster and gene tested
        mutate(ubc_subcluster = clust, .before = 1) %>%
        rownames_to_column(var = "gene")

      # return differential expression results
      return(de_genes)
    }
  ) %>%
    # combine all results into one dataframe
    do.call(what = rbind, args = .)

  return(de_genes)
}

prep_and_make_volcano <- function(
  de_genes,
  logfc_column,
  pval_column,
  genes_to_label,
  plot_title = NULL
) {
  # get direction
  de_genes <- mutate(
    de_genes,
    direction = case_when(
      rlang::parse_expr(
        paste0(logfc_column, " < 0 & ", pval_column, " < 0.05 ~ 'down'")
      ),
      rlang::parse_expr(
        paste0(logfc_column, " > 0 & ", pval_column, " < 0.05 ~ 'up'")
      ),
      .default = "n.s."
    ),
    gene_label = case_when(gene %in% genes_to_label ~ gene, .default = "")
  )

  # show number of up/downregulated genes on volcano plot
  caption <- sprintf(
    "n = %s; down = %s; up = %s",
    nrow(de_genes),
    sum(de_genes$direction == "down"),
    sum(de_genes$direction == "up")
  )

  # make volcano plot
  .plt <- make_volcano(
    data = de_genes,
    log_fc = .data[[logfc_column]],
    log_pval = -log10(.data[[pval_column]]),
    direction = direction,
    gene_label = gene_label,
    label_num_genes = FALSE,
    seed = 256
  ) + 
    labs(
      title = plot_title,
      caption = caption,
      x = "log2(FC)",
      y = "-log10(p-value)"
    ) + 
    scale_colour_manual(
      values = c(down = "red", up = "blue", `n.s.` = "black")
    ) + 
    theme_classic() + 
    theme(
      axis.text = element_text(colour = "black"),
      legend.position = "none"
    )

  return(.plt)
}

# ------------------------------------------------------------------------------
# load Seurat object

# human UBCs with clustering by Quang
srat <- readRDS(arg_list$srat_file)
srat$dataset_name <- str_remove(srat$dataset_name, "_full_cerebellum_human")

# ------------------------------------------------------------------------------
# prepare for differential gene expression analysis

# split Seurat object by dataset
srat_list <- SplitObject(srat, split.by = "dataset_name")

# which genes are expressed in >= 1% of cells in both datasets?
genes_to_keep <- map(
  .x = names(srat_list),
  .f = \(x) {
    gene_is_expressed <- GetAssayData(object = srat_list[[x]], slot = "counts", assay = "RNA") > 0
    keep_gene <- Matrix::rowSums(gene_is_expressed) >= (0.01 * ncol(srat_list[[x]]))
    message(sprintf("***Keeping %s genes from the %s dataset***", sum(keep_gene), x))

    # return vector of genes to keep
    return(names(keep_gene)[keep_gene])
  }
) %>%
  setNames(nm = names(srat_list))
genes_to_keep$intersect <- purrr::reduce(.x = genes_to_keep, .f = intersect)
saveRDS(genes_to_keep, file = file.path(arg_list$out_dir, "genes_to_keep.rds"))

message(sprintf(
  "***Keeping %s genes expressed in >= 1%% of cells in both datasets***",
  length(genes_to_keep$intersect)
))
srat_list <- map(
  .x = srat_list,
  .f = \(x) {
    # subset for expressed genes
    x <- subset(x, features = genes_to_keep$intersect) %>%
      # re-run SCTransform
      SCTransform(vars.to.regress = "CC.Difference", variable.features.n = 3000)
  }
)

# merge datasets again
srat <- merge(x = srat_list[[1]], y = srat_list[-1])

# prep for differential gene expression
srat <- PrepSCTFindMarkers(srat)

# ------------------------------------------------------------------------------
# differential gene expression

# set levels for clusters and heatmap afterwards
srat$ubc_subcluster <- factor(paste0("UBC_", srat$SCT_snn_res.0.5))
srat$age <- factor(srat$age, levels = str_sort(unique(srat$age), numeric = TRUE))

# set Idents to UBC subclusters
Idents(srat) <- srat$ubc_subcluster

if (arg_list$conserved) {
  # run FindConservedMarkers
  de_genes <- run_all_de(
    srat = srat,
    clust = levels(srat$ubc_subcluster)
  )

  # rank top genes by log2FC
  feature_ranking <- list()
  for (clust in levels(srat$ubc_subcluster)) {
    feature_ranking[[clust]] <- de_genes %>%
      # get genes that are upregulated in both datasets
      filter(
        ubc_subcluster == clust,
        Aldinger_p_val < 0.05,
        Sepp_p_val < 0.05,
        Aldinger_avg_log2FC > 0,
        Sepp_avg_log2FC > 0
      ) %>%
      # get each feature's rank in each dataset
      mutate(
        aldinger_rank = min_rank(desc(Aldinger_avg_log2FC)),
        sepp_rank = min_rank(desc(Sepp_avg_log2FC)),
        rank = aldinger_rank + sepp_rank
      ) %>%
      arrange(rank) %>%
      select(gene, aldinger_rank, sepp_rank, rank)
  }
} else {
  # run FindMarkers
  de_genes <- run_all_de(
    srat = srat,
    clust = levels(srat$ubc_subcluster)
  )

  # rank top genes by log2FC
  feature_ranking <- list()
  for (clust in levels(srat$ubc_subcluster)) {
    feature_ranking[[clust]] <- de_genes %>%
      filter(ubc_subcluster == clust, p_val < 0.05, avg_log2FC > 0) %>%
      mutate(rank = min_rank(desc(avg_log2FC))) %>%
      arrange(rank) %>%
      select(gene, rank)
  }
}

# save differential gene expression results
write_tsv(
  x = de_genes,
  file = file.path(arg_list$out_dir, "de_genes.tsv")
)

# save ranking of upregulated genes
saveRDS(
  feature_ranking,
  file = file.path(arg_list$out_dir, "cluster_marker_ranking.rds")
)

# ------------------------------------------------------------------------------
# volcano plots

walk(
  .x = levels(srat$ubc_subcluster),
  .f = \(clust) {
    message(sprintf("***Generating volcano plot for %s***", clust))

    # subset for cluster
    clust_de <- de_genes[de_genes$ubc_subcluster == clust, ]

    # make plot for each dataset
    if (arg_list$conserved) {
      .plt <- map(
        .x = c("Aldinger", "Sepp"),
        .f = \(dataset) {
          # columns to use
          logfc_column <- paste0(str_to_title(dataset), "_avg_log2FC")
          pval_column <- paste0(str_to_title(dataset), "_p_val")

          # get top genes (by log2FC) to label in each dataset
          top_genes <- union(
            x = slice_max(clust_de, Aldinger_avg_log2FC, n = 10) %>% pull("gene"),
            y = slice_max(clust_de, Sepp_avg_log2FC, n = 10) %>% pull("gene")
          )

          # process DE results and make volcano
          .plt <- prep_and_make_volcano(
            de_genes = clust_de,
            logfc_column = logfc_column,
            pval_column = pval_column,
            genes_to_label = top_genes,
            plot_title = paste0(clust, " (", dataset, ")")
          )
          return(.plt)
        }
      )

      # combine Aldinger/Sepp plots
      .plt <- wrap_plots(.plt, nrow = 1)

      # plot params
      width = 12
      height = 6
    } else {
      # get top genes by log2FC
      top_genes <- slice_max(clust_de, avg_log2FC, n = 10) %>% pull("gene")

      # process DE results and make volcano
      .plt <- prep_and_make_volcano(
        de_genes = clust_de,
        logfc_column = "avg_log2FC",
        pval_column = "p_val",
        genes_to_label = top_genes,
        plot_title = clust
      )

      # plot params
      width = 6
      height = 6
    }

    # save plots
    ggsave(
      filename = paste0(clust, "_volcano.png"),
      plot = .plt,
      path = arg_list$out_dir,
      width = width,
      height = height,
      units = "in",
      dpi = 600
    )
  }
)

# ------------------------------------------------------------------------------
# heatmap

# get average age of each cluster to show them youngest --> oldest
ubc_by_age <- srat[[c("age", "ubc_subcluster")]] %>%
  # convert age to number
  mutate(age = as.numeric(str_remove(string = age, pattern = " PCW"))) %>%
  # get mean age for each ubc_subcluster
  group_by(ubc_subcluster) %>%
  summarise(avg_age = mean(age)) %>%
  # sort by age
  arrange(avg_age) %>%
  # extract ubc_subcluster
  pull(ubc_subcluster) %>%
  as.character()

# create top annotation rows
ages <- srat[[]] %>%
  select(ubc_subcluster, age) %>%
  table() %>%
  prop.table(margin = 1) %>%
  .[ubc_by_age, ]
sex <- srat[[]] %>%
  select(ubc_subcluster, sex) %>%
  table() %>%
  prop.table(margin = 1) %>%
  .[ubc_by_age, ]
dset <- srat[[]] %>%
  select(ubc_subcluster, dataset_name) %>%
  mutate(dataset_name = str_remove(dataset_name, "_full_cerebellum_human")) %>%
  table() %>%
  prop.table(margin = 1) %>%
  .[ubc_by_age, ]
pal1 <- scales::pal_gradient_n(
  scales::pal_brewer(palette = "YlOrRd")(3)
)(
  seq(0, 1, length.out = ncol(ages))
)
pal2 <- c("pink","blue") # sex
pal3 <- RColorBrewer::brewer.pal(3, "Dark2") # dataset

# save top annotation rows
walk2(
  .x = list(ages, sex, dset),
  .y = c("age", "sex", "dataset"),
  .f = \(x, y) {
    saveRDS(x, file.path(arg_list$out_dir, paste0("prop_", y, ".rds")))
  }
)

topannot <- HeatmapAnnotation(
  age = anno_barplot(
    ages, 
    gp = gpar(fill = pal1, border="white", lty="solid"),
    bar_width = 1, 
    height = unit(2, "in"),
    axis=FALSE
  ),
  sex = anno_barplot(
    sex,
    gp = gpar(fill = pal2, border="white"),
    bar_width = 1,
    height = unit(1,"in"),
    axis=FALSE
  ),
  dataset = anno_barplot(
    dset,
    gp = gpar(fill = pal3, border="white"),
    bar_width = 1,
    height = unit(1,"in"),
    axis=FALSE
  )
)
lgd_list <- list(
  Legend(labels = colnames(ages), title = "age", legend_gp = gpar(fill = pal1)),
  Legend(labels = colnames(sex), title = "sex", legend_gp = gpar(fill = pal2)),
  Legend(labels = colnames(dset), title = "dataset", legend_gp = gpar(fill = pal3))
)

# top markers to show on heatmap
unique_markers <- map(
  .x = ubc_by_age,
  .f = \(clust) {
    x <- feature_ranking[[clust]]
    nr <- min(nrow(x), 100)
    x <- pull(x, gene) %>%
      head(nr)
    return(x)
  }
) %>%
  unlist() %>%
  unique()

# get data slot, then manually scale the results
xpr <- AverageExpression(
  srat,
  assays = "SCT",
  features = unique_markers,
  group.by = "ubc_subcluster",
  slot = "data"
)[[1]][, ubc_by_age]
xpr <- t(scale(t(xpr)))

# heatmap colours
col_fun <- circlize::colorRamp2(
  breaks = c(-2,0,2),
  colors = c("#FF00FF","#000000", "#FFFF00")
)

# label top genes
anno_genes <- map(
  .x = feature_ranking,
  .f = \(df) {slice_min(df, order_by = rank, n = 5)}
) %>%
  unlist() %>%
  unname() %>%
  c("OTX2", "CBFA2T2", "SOX4", "SOX11", "LMX1A", "EOMES")
ha <- rowAnnotation(
  foo = anno_mark(at = match(anno_genes, rownames(xpr)), labels = anno_genes)
)

hm <- Heatmap(
  xpr,
  name="scaled\nexpression",
  cluster_rows=FALSE,
  cluster_columns=FALSE,
  col = col_fun,
  show_row_names = FALSE,
  top_annotation = topannot,
  right_annotation = ha
)

png(
  filename = file.path(arg_list$out_dir, "avg_expr_heatmap.png"),
  width = 12,
  height = 16,
  units = "in",
  res = 300
)
draw(hm, annotation_legend_list = lgd_list, merge_legend = TRUE)
dev.off()

# ------------------------------------------------------------------------------

message("\n***SESSION INFO***\n")
print(sessionInfo())

