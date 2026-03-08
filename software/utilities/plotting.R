#' Load plotting functions from other files.
#'
#' @param d Directory containing `files`.
#' @param files Vector of file names to load (`source`).
#'
#' @return
#'
load_viz_functions <- function(
  d = "/home/rstudio/isilon/private/icheong/CBL_scRNAseq/software/utilities/",
  files = c("cluster_barplot.R", "plot_venn_diagrams.R")
) {
  full_paths <- file.path(d, files)
  walk(
    .x = full_paths,
    .f = \(x) {
      message(sprintf("Loading functions from %s", x))
      source(x)
    }
  )
}


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


#' Make DimPlots highlighting each group in a metadata column.
#' 
#' @param srat Seurat object.
#' @param highlight_by The metadata column to generate plots for. For example,
#'   to generate plots highlighting cells from each cluster, use `highlight_by =
#'   "seurat_clusters"`.
#' @param highlight_cols Colours to highlight the cells as. Can be a single
#'   colour or vector of colours. If the vector of colours is less than the
#'   number of groups in `highlight_by`, then the colours will be repeated.
#' @param ... Additional parameters passed to `DimPlot`.
#' 
#' @return A list of ggplot objects.
#' 
highlight_DimPlot <- function(
  srat,
  highlight_by,
  highlight_cols = NULL,
  ...
) {
  # get vector of cell clusters to highlight
  if (is.factor(srat[[]][, highlight_by])) {
    clusters <- levels(srat[[]][, highlight_by])
  } else {
    clusters <- str_sort(unique(srat[[]][, highlight_by]))
  }

  # set colours
  if (is.null(highlight_cols)) {
    highlight_cols <- scales::pal_hue()(length(clusters))
  } else if (length(highlight_cols) < length(clusters)) {
    highlight_cols <- rep(highlight_cols, length.out = length(clusters))
  }

  # get list of plots highlighting each cell cluster
  .plt_list <- map2(
    .x = clusters,
    .y = highlight_cols,
    .f = \(x, y) {
      .plt <- highlight_SingleDimPlot(
        srat = srat,
        highlight_by = highlight_by,
        highlight_group = x,
        cols.highlight = y,
        ...
      ) + 
        NoLegend() + 
        NoAxes() + 
        labs(title = x)

      return(.plt)
    }
  )

  return(.plt_list)
}


#' Make a single DimPlot highlighting a group in a metadata column. Uses
#' Seurat's `cells.highlight` from `DimPlot`.
#' 
#' @param srat Seurat object.
#' @param highlight_by The metadata column to generate plots for. For example,
#'   to generate plots highlighting cells from each cluster, use `highlight_by =
#'   "seurat_clusters"`.
#' @param highlight_group The group or vector of groups from the metadata column
#'   to highlight.
#' @param ... Additional parameters passed to `DimPlot`.
#' 
#' @return A single ggplot object.
#' 
highlight_SingleDimPlot <- function(
  srat,
  highlight_by,
  highlight_group,
  ...
) {
  # get cells to highlight
  cells_highlight <- filter(
    srat[[]],
    .data[[highlight_by]] %in% highlight_group
  ) %>%
    rownames()

  # make plot
  .plt <- do.call(
    what = DimPlot,
    args = list(
      object = srat,
      cells.highlight = cells_highlight,
      ...
    )
  )

  return(.plt)
}



load_viz_functions()
