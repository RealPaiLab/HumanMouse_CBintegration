# ==============================================================================
# Run differential gene expression on integrated human UBC subclusters.
# Copied over on 250417. from commit 39938e6
# ==============================================================================

library(argparse)
library(tidyverse)
library(patchwork)
library(Seurat)
library(ComplexHeatmap)

#srat_file <- "/home/rstudio/isilon/public/HumanMouseUBC/data/UBC.Harmony.RDS"
srat_file <- "/home/rstudio/isilon/private/projects/HumanMouseUBC/UBCclusters/fromQuang/UBC.Harmony.RDS"
test_use <- "MAST"
latent_vars <- c("sex","dataset_name")
out_dir <- "/home/rstudio/isilon/private/projects/HumanMouseUBC/integrated_human_ubc/diffExpr"

useDEfile <- sprintf("%s/250424/de_genes.tsv", out_dir)
useFRfile <- sprintf("%s/250424/cluster_marker_ranking.rds", out_dir)

colourFile <- "../UBCcolours.txt"
clrs <- read.delim(colourFile,header=TRUE,sep="\t")

heatmap_logFC_thresh <- 1.5 # logFC threshold for heatmap
heatmap_pval_thresh <- 0.05 # adjusted p-value threshold for heatmap

dt <- format(Sys.Date(),"%y%m%d")
out_dir <- sprintf("%s/%s",out_dir,dt)

if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = FALSE)
}

logFile <- sprintf("%s/log.txt", out_dir)
sink(logFile,split=TRUE)
tryCatch({

# ------------------------------------------------------------------------------
# functions

#' Make a volcano plot. Assumes that the values have already been
#' log-transformed.
#'
#' @param data Data to be plotted.
#' @param log_fc Expression (not string) containing column name of the log
#'   fold-changes in the data.
#' @param log_pval Expression (not string) containing column name of the log
#'   p-values in the data.
#' @param direction Expression (not string) containing column name of the
#'   expression direction (e.g. up or down).
#' @param gene_labels Expression (not string) containing column name in data
#'   with the gene names to add to the plot.
#' @param log_fc_thresh Add lines marking a log fold-change cutoff. For example,
#'   to show an absolute logFC greater than 1, use `c(1, -1)`.
#' @param log_pval_thresh Same as above but for the log p-value. Note that this
#'   should be a log-transformed value, e.g. `-log10(0.05)`.
#' @param label_num_genes Add a line showing the total number of genes, the
#'   number of upregulated genes, and the number of downregulated genes. If
#'   `TRUE` (default), `log_pval_threshold` cannot be `NULL`.
#' @param seed Set seed for the gene labels (passed to
#'   `ggrepel::geom_text_repel`). Defaults to NA.
#'
#' @return A `ggplot` object.
#'
make_volcano <- function(
  data,
  log_fc,
  log_pval,
  direction,
  gene_labels = NULL,
  log_fc_thresh = NULL,
  log_pval_thresh = NULL,
  label_num_genes = TRUE,
  seed = NA
) {
  plt <- ggplot(
    data = data,
    mapping = aes(
      x = {{ log_fc }},
      y = {{ log_pval }},
      colour = {{ direction }}
    )
  ) + 
    geom_point(size = 1) + 
    labs(x = "log FC", y = "-log p-value") + 
    theme_minimal() + 
    theme(axis.text = element_text(colour = "black"))
  
  # add labels for gene names
  if (!is.null(eval(substitute(gene_labels), data))) {
    plt <- plt + 
      ggrepel::geom_text_repel(
        aes(label = {{ gene_labels }}),
        max.overlaps = Inf,
        show.legend = FALSE,
        seed = seed
      )
  }
  
  # add vertical line for logFC
  if (!is.null(log_fc_thresh)) {
    plt <- plt + 
      geom_vline(
        xintercept = log_fc_thresh,
        colour = "red",
        linewidth = 0.25,
        linetype = "dashed"
      )
  }
  
  # add horizontal line for p-value
  if (!is.null(log_pval_thresh)) {
    plt <- plt + 
      geom_hline(
        yintercept = log_pval_thresh,
        colour = "red",
        linewidth = 0.25,
        linetype = "dashed"
      )
  }
  
  # add number of (differentially expressed) genes to plot
  if (label_num_genes) {
    df <- mutate(
      data,
      log_fc = {{ log_fc }},
      log_pval = {{ log_pval }},
      .keep = "none"
    )
    
    n_up <- sum(df$log_fc > 0 & df$log_pval > log_pval_thresh)
    n_down <- sum(df$log_fc < 0 & df$log_pval > log_pval_thresh)
    
    label <- paste0(
      "n = ", nrow(data),
      "; down = ", n_down,
      "; up = ", n_up
    )
    plt <- plt + labs(subtitle = label)
  }
  
  return(plt)
}

#' Add column named of gene labels for downstream volcano plot.
#'
#' @param dge_res A `TopTags` object or a `data.frame`.
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

# ' Run differential gene expression for all clusters.
run_all_de <- function(
  srat,
  clusters,
  test_use = test_use
) {
  # run differential gene expression (each cluster vs. all other clusters)
  de_genes <- map(
    .x = clusters,
    .f = \(clust) {
      message(sprintf("***%s, n = %s cells***", clust, sum(Idents(srat) == clust)))

      # differential expression
      de_genes <- FindMarkers(
        object = srat,
        ident.1 = clust,
        grouping.var = "dataset_name",
        assay = "SCT",
        logfc.threshold = 0,
        min.pct = 0.01,
        test.use = "MAST",
        latent.vars = c("sex","dataset_name")
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

# ################################################################################# Work starts here.

cat("Reading human UBC cluster file\n")
t0 <- Sys.time()
srat <- readRDS(srat_file)
print(Sys.time() - t0)
# print number of cells and genes
cat(sprintf("Number of cells: %s\n", ncol(srat)))
cat(sprintf("Number of genes: %s\n", nrow(srat)))
srat$dataset_name <- str_remove(srat$dataset_name, "_full_cerebellum_human")

# remove samples with age "9 PCW" or "10 PCW"
srat <- subset(srat, subset = !age %in% c("9 PCW", "10 PCW"))

# split Seurat object by dataset
srat_list <- SplitObject(srat, split.by = "dataset_name")

cat("Identifying genes expressed in >= 1%% of cells in both datasets\n")
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

saveRDS(genes_to_keep, file = file.path(out_dir, "genes_to_keep.rds"))

message(sprintf(
  "***Keeping %s genes expressed in >= 1%% of cells in both datasets***",
  length(genes_to_keep$intersect)
))

cat("Re-running SCTransform on each dataset using intersecting genes\n")
t0 <- Sys.time()
srat_list <- map(
  .x = srat_list,
  .f = \(x) {
    # subset for expressed genes
    x <- subset(x, features = genes_to_keep$intersect) %>%
      # re-run SCTransform
      SCTransform(
        vars.to.regress = "CC.Difference", 
        variable.features.n = 3000)
  }
)
print(Sys.time() - t0)

# merge datasets again
srat <- merge(x = srat_list[[1]], y = srat_list[-1])

# prep for differential gene expression
cat("Running PrepSCTFindMarkers\n")
t0 <- Sys.time()
srat <- PrepSCTFindMarkers(srat)
print(Sys.time() - t0)

# set levels for clusters and heatmap afterwards
srat$ubc_subcluster <- factor(paste0("UBC_", srat$SCT_snn_res.0.5))
srat$age <- factor(srat$age, 
    levels = str_sort(unique(srat$age), numeric = TRUE))

# set Idents to UBC subclusters
Idents(srat) <- srat$ubc_subcluster

if (file.exists(useDEfile)) {
  cat("DE file exists. loading\n")
  de_genes <- read.delim(file = useDEfile,sep="\t",h=T,as.is=T)
  feature_ranking <- readRDS(file = useFRfile)
} else {
  # run FindMarkers
  cat("Running FindMarkers\n")
  t0 <- Sys.time()
    de_genes <- run_all_de(
      srat = srat,
      clust = levels(srat$ubc_subcluster)
  )
  print(Sys.time() - t0)

    # rank top genes by log2FC
    feature_ranking <- list()
    for (clust in levels(srat$ubc_subcluster)) {
      feature_ranking[[clust]] <- de_genes %>%
        filter(ubc_subcluster == clust, 
          p_val < heatmap_pval_thresh, 
          avg_log2FC > heatmap_logFC_thresh
        ) %>%
        mutate(rank = min_rank(desc(avg_log2FC))) %>%
        arrange(rank) %>%
        select(gene, rank)
  }

  # save differential gene expression results
  write_tsv(
    x = de_genes,
    file = file.path(out_dir, "de_genes.tsv")
  )

  cat("save ranking of upregulated genes\n")
  saveRDS(
    feature_ranking,
    file = file.path(out_dir, "cluster_marker_ranking.rds")
  )
}

# ------------------------------------------------------------------------------
# volcano plots
walk(
  .x = levels(srat$ubc_subcluster),
  .f = \(clust) {
    cat(sprintf("***Generating volcano plot for %s***", clust))

    # subset for cluster
    clust_de <- de_genes[de_genes$ubc_subcluster == clust, ]

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

    # save plots
    ggsave(
      filename = paste0(clust, "_volcano.png"),
      plot = .plt,
      path = out_dir,
      width = width,
      height = height,
      units = "in",
      dpi = 600
    )
  }
)
cat("Done saving volcano plots.\n")

# ------------------------------------------------------------------------------
# heatmap

cat("Removing cluster 6 ******\n")
srat <- subset(srat, idents = c("UBC_0", "UBC_1", "UBC_2", "UBC_3", "UBC_4", "UBC_5", "UBC_7", "UBC_8"))

cat("* Generating heatmap\n")
cat("\t Creating top annotation rows\n")
browser()
# get average age of each cluster to show them youngest --> oldest
###ubc_by_age <- srat[[c("age", "ubc_subcluster")]] %>%
###  # convert age to number
###  mutate(age = as.numeric(
###    str_remove(string = age, pattern = " PCW"))) %>%
###  # get mean age for each ubc_subcluster
###  group_by(ubc_subcluster) %>%
###  summarise(avg_age = mean(age)) %>%
###  # sort by age
###  arrange(avg_age) %>%
###  # extract ubc_subcluster
###  pull(ubc_subcluster) %>%
###  as.character()

clrs$Cluster <- paste("UBC_", clrs$Cluster, sep = "")
clrs <- clrs[-which(clrs$Cluster %in% "UBC_6"),]
ubc_by_age <- clrs$Cluster
srat$ubc_subcluster <- factor(srat$ubc_subcluster,levels=ubc_by_age)


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
    saveRDS(x, file.path(out_dir, paste0("prop_", y, ".rds")))
  }
)

lgd_format <- gpar(fontsize=28)
lgd_ttl <- gpar(fontsize=28, fontface = "bold")

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
    annotation_name_gp = gpar(fontsize=28, fontface = "bold"),
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
  Legend(labels = colnames(ages), title = "age", 
    legend_gp = gpar(fill = pal1), labels_gp = lgd_format, title_gp = lgd_ttl),
  Legend(labels = colnames(sex), title = "sex", 
    legend_gp = gpar(fill = pal2), labels_gp = lgd_format, title_gp = lgd_ttl),
  Legend(labels = colnames(dset), title = "dataset", 
    legend_gp = gpar(fill = pal3), labels_gp = lgd_format, title_gp = lgd_ttl)
)
cat("\t\t Done creating top annotation rows\n")

cat("\t Creating heatmap data\n")
# top markers to show on heatmap
cat("* Fetching unique markers for heatmap plotting\n")
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

cat(sprintf("%i unique markers\n", length(unique_markers)))

browser()

cat("Getting average expression for cluster markers\n")
# get data slot, then manually scale the results
xpr <- AverageExpression(
  srat,
  assays = "SCT",
  features = unique_markers,
  group.by = "ubc_subcluster",
  slot = "data"
)[[1]]
colnames(xpr) <- sub("-","_", colnames(xpr))
xpr <- xpr[, ubc_by_age]
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
  foo = anno_mark(at = match(anno_genes, rownames(xpr)), labels = anno_genes,
    labels_gp=gpar(fontsize = 24, col = "black"))
)

# order by cluster
hm <- Heatmap(
  xpr,
  name="scaled\nexpression",
  cluster_rows=FALSE,
  cluster_columns=FALSE,
  col = col_fun,
  show_row_names = FALSE,
  top_annotation = topannot,
  right_annotation = ha,
  column_names_rot = 45,
  column_names_gp = gpar(fontsize = 20),
  heatmap_legend_param = list(labels_gp = lgd_format, title_gp = lgd_ttl)
)

pdf(
  file = file.path(out_dir, "avg_expr_heatmap.pdf"),
  width = 14,
  height = 20
)
draw(hm, annotation_legend_list = lgd_list, merge_legend = TRUE)
dev.off()

# ------------------------------------------------------------------------------

message("\n***SESSION INFO***\n")
print(sessionInfo())

},error=function(ex){
    print(ex)
},
finally={
    sink()
}
)