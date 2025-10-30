#' Run differential expression and save genes to CSV file.
#'
#' @param cds Monocle `cell_data_set` object
#' @param neighbor_graph Either "knn" or "principal_graph", passed to
#'   `graph_test`
#' @param cores Number of cores to use
#' @param deg_file Location of output file containing the differential gene
#'   expression results (full path required)
#'
#' @return A dataframe of differentially expressed genes (output of
#'   `graph_test`)
#' 
run_graph_test <- function(
  cds,
  neighbor_graph,
  cores,
  deg_file = NULL
) {
  message("***Running differential expression with `graph_test`")
  gt_res <- graph_test(
    cds,
    neighbor_graph = neighbor_graph,
    cores = cores
  )
  
  if (!is.null(deg_file)) {
    message(sprintf("***Saving differentially expressed genes to %s", deg_file))
    write_csv(
      gt_res %>% rownames_to_column(),
      file = deg_file
    )
  }
  
  return(gt_res)
}


#' Find gene modules and save them to CSV file.
#'
#' @param cds Monocle `cell_data_set` object
#' @param resolution Passed to `find_gene_modules`
#' @param seed Set random seed
#' @param verbose Print verbose output
#' @param module_file Location of output file containing the genes associated
#'   with each module (full path required)
#'
#' @return A dataframe of genes and their associated modules (output of
#'   `find_gene_modules`)
#' 
get_modules <- function(
  cds,
  resolution,
  seed = 42,
  verbose = TRUE,
  module_file = NULL
) {
  message("***Finding gene modules with `find_gene_modules`")
  modules_df <- find_gene_modules(
    cds = cds,
    resolution = resolution,
    random_seed = seed,
    verbose = verbose
  )
  
  if (!is.null(module_file)) {
    message(sprintf("***Saving gene modules to %s", module_file))
    write_csv(modules_df, file = module_file)
  }
  
  return(modules_df)
}


#' Get aggregated expression of gene modules and save to CSV file
#'
#' @param cds Monocle `cell_data_set` object
#' @param modules_df Dataframe of genes in the first column and their module in
#'   the second column
#' @param cell_annot_df Dataframe of cell IDs in the first column and their
#'   annotation in the second column (e.g., cell type, age, etc.)
#' @param out_dir Directory to save the aggregated gene module expression matrix
#'
#' @return Matrix of gene module expression
#' 
get_module_expr <- function(
  cds,
  modules_df,
  cell_annot_df,
  out_dir
) {
  module_expr <- aggregate_gene_expression(
    cds = cds,
    gene_group_df = modules_df,
    cell_group_df = cell_annot_df
  )
  
  module_expr_file <- file.path(
    out_dir,
    paste0("module_expr_", colnames(cell_annot_df)[2], ".csv")
  )
  message(sprintf("***Saving gene module expression matrix to %s", module_expr_file))
  write.csv(module_expr, file = module_expr_file)
  
  # base::rownames(module_expr) <- str_c("Module ", base::rownames(module_expr))
  return(module_expr)
}


#' Make heatmap of aggregated gene module expression and save it.
#'
#' @param module_expr Matrix of gene module expression
#' @param cell_annot Name of cell annotation column (to be included in the
#'   filename of the heatmap)
#' @param out_dir Directory to save heatmap plot to
#'
#' @return `pheatmap` plot
#' 
plot_gene_modules <- function(
  module_expr,
  cell_annot,
  out_dir
) {
  message(sprintf("***Plotting gene module expression heatmap for %s", cell_annot))
  hm <- pheatmap::pheatmap(
    mat = module_expr,
    color = scales::viridis_pal()(100),
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    scale = "column",
    clustering_method = "ward.D2",
    silent = TRUE
  )
  
  hm_file <- paste0("heatmap_", cell_annot, ".png")
  message(sprintf("***Saving gene module expression heatmap to %s", hm_file))
  
  ggsave(
    hm_file,
    plot = hm,
    path = out_dir,
    width = 0.25 * dim(module_expr)[2] + 1,
    height = 0.25 * dim(module_expr)[1] + 1,
    units = "in",
    dpi = 600
  )
  
}


#' Plot a UMAP from a `cell_data_set` object and highlight specific cells.
#' 
#' Similar to `cells.highlight` in Seurat's `DimPlot` function, but this is
#'   meant for Monocle 3 which doesn't have an equivalent function (at least not
#'   one that I know of).
#'
#' @param cds A `cell_data_set` object from Monocle 3.
#' @param cells A **named** list of cells to highlight. The names of the list
#'   will be the label on the plot.
#' @param cell_size The size of the cells on the plot.
#' @param highlight_cols A vector of colours to highlight the cells; must be the
#'   same length as `cells`. If `NULL` (default), the cells will be highlighted 
#'   in with the default `ggplot2` colours.
#' @param na_col The colour of the non-highlighted cells.
#'
#' @return A `ggplot2` object.
#'
plot_highlight <- function(
  cds,
  cells,
  cell_size = 0.1,
  highlight_cols = NULL,
  na_col = "grey"
) {
  # extract the UMAP embeddings
  plot_df <- reducedDim(x = cds, type = "UMAP") %>%
    as.data.frame() %>%
    dplyr::rename(UMAP_1 = V1, UMAP_2 = V2)
  
  # convert list of cells to dataframe of cells to highlight
  hl <- data.frame(
    highlight = rep(x = names(cells), times = purrr::map(cells, length)),
    row.names = unlist(cells, use.names = FALSE)
  ) %>%
    mutate(highlight = factor(highlight, levels = names(cells)))

  # add to plotting
  plot_df <- merge(x = plot_df, y = hl, by = "row.names", all.x = TRUE) %>%
    dplyr::rename(cell_ids = Row.names) %>%
    # sort rows so NA values are plotted first
    arrange(!is.na(highlight))

  # make plot
  .plt <- ggplot(plot_df, aes(x = UMAP_1, y = UMAP_2, colour = highlight)) + 
    geom_point(size = cell_size) + 
    theme_classic() + 
    theme(
      axis.text = element_text(colour = "black"),
      axis.ticks = element_line(colour = "black"),
      legend.title = element_blank()
    )

  # set the colours of the highlighted cells
  if (is.null(highlight_cols)) {
    # default `ggplot2` colours if no custom colours are provided
    .plt <- .plt + 
      scale_colour_hue(
        na.value = na_col,
        guide = guide_legend(override.aes = list(size = 2))
      )
  } else {
    # custom colours
    .plt <- .plt + 
      scale_colour_manual(
        values = highlight_cols,
        na.value = na_col,
        guide = guide_legend(override.aes = list(size = 2))
      )
  }

  return(.plt)
}


#' Extract and save pseudotime values from a `cell_data_set` object
#'
#' @param cds A `cell_data_set` object from Monocle 3.
#' @param file File to write to.
#'
save_pseudotime_to_csv <- function(
  cds,
  file
) {
  # extract pseudotime values
  ptime <- data.frame(pseudotime = pseudotime(cds)) %>%
    rownames_to_column("cell_id")

  # save to CSV
  message(sprintf("Saving `pseudotime` values to %s", file))
  write_csv(
    x = ptime,
    file = file
  )
}
