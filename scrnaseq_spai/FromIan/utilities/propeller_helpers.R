#' Get and format the cell type proportions in `propeller` for plotting with
#' `plot_props_per_sample`.
#'
#' @param clusters A factor specifying the cluster or cell type for each cell
#'   (passed to `speckle::getTransformedProps`).
#' @param sample A factor specifying the biological replicate for each cell
#'   (passed to `speckle::getTransformedProps`).
#' @param sample_meta_df A dataframe with sample metadata to be added to the output (e.g., dataest, species, age, etc.). Must contain a column named "sample". Optional argument (defaults to NULL).
#'
#' @return A dataframe containing the output of `propeller::getTransformedProps`
#'   with the proportion of each cluster/cell type in each sample and
#'   (optionally) additional sample metadata if provided.
#'
get_props_per_sample <- function(
  clusters,
  sample,
  sample_meta_df = NULL
) {
  props <- getTransformedProps(
    clusters = clusters,
    sample = sample,
    transform = "logit"
  ) %>%
    pluck("Proportions") %>% # get this list item
    as.data.frame() %>%
    rename(proportion = Freq)

  if (!is.null(sample_meta_df)) {
    props <- merge(x = props, y = sample_meta_df, by = "sample", all.x = TRUE)
  }

  return(props)
}


#' Make bar plot showing the proportion of cells in each cluster for each
#' sample.
#'
#' @param props The output of `get_props_per_sample` (a dataframe with the
#'   proportion of cells in each cluster/cell type for each sample).
#' @param facet_cols Name of variable for which to create a faceting column.
#'
#' @return A ggplot2 object.
#'
plot_props_per_sample <- function(
  props,
  facet_cols
) {
  # make the plot
  .plt <- ggplot(
    props,
    aes(x = sample, y = proportion, fill = clusters)
  ) + 
    geom_col(position = "stack") + 
    facet_grid(cols = vars({{facet_cols}}), scales = "free", space = "free") + 
    scale_y_continuous(expand = expansion(0, 0.05)) + 
    scale_fill_manual(
      values = pals::cols25(n = length(unique(props[["clusters"]])))
    ) + 
    theme_minimal() + 
    theme(
      axis.text.x = element_text(hjust = 1, angle = 60)
    )

  return(.plt)
}


#' Get the FDRs for each cluster after running propeller.
#'
#' @param propeller_results The output of `speckle::propeller`.
#'
#' @return A dataframe containing with the cluster and its FDR.
#'
get_fdrs_for_plotting <- function(
  propeller_results
) {
  fdrs <- propeller_results[, c("BaselineProp.clusters", "FDR")] %>%
    rename(clusters = BaselineProp.clusters) %>%
    arrange(clusters) %>%
    remove_rownames()
  return(fdrs)
}


#' Make boxplot showing proportion of cells in each cluster.
#'
#' @param props The output of `get_props_per_sample` (a dataframe with the
#'   proportion of cells in each cluster/cell type for each sample).
#' @param facet_var Name of variable (passed to `facet_wrap`).
#' @param facet_ncol Number of columns (passed to `facet_wrap`).
#' @param fdrs The output of `get_fdrs_for_plotting` (a data frame with the FDRs
#'   for each cluster).
#' @param fdr_y The y-value at which to plot the FDR values. If `NULL`
#'   (default), the y-values will be calculated automatically for each facet.
#' @param fdr_label_size Size of the FDR labels (in "points"). Passed to
#'   `geom_text`. Defaults to `11`.
#'
#' @return A `ggplot2` object.
#'
boxplot_props_per_cluster <- function(
  props,
  facet_var,
  facet_ncol = 4,
  fdrs = NULL,
  fdr_y = NULL,
  fdr_label_size = 11
) {
  # 
  .plt <- ggplot(props, aes(x = species, y = proportion)) + 
    geom_boxplot(width = 0.6, outlier.shape = NA) + 
    geom_jitter(aes(shape = dataset_name), width = 0.2)
  
  if (!is.null(fdrs)) {
    # where to show the FDR value (how high on the plot)
    if (is.null(fdr_y)) {
      fdrs <- props %>%
        group_by(clusters) %>%
        summarise(fdr_y = max(proportion) * 1.05) %>%
        full_join(
          x = fdrs,
          y = .,
          by = join_by(clusters)
        )
    }

    # show FDR value on plot
    .plt <- .plt + 
      geom_text(
        data = fdrs,
        mapping = aes(x = "mouse", y = fdr_y, label = sprintf("adj. p = %.3f", FDR)),
        nudge_x = 0.5,
        size = fdr_label_size,
        size.unit = "pt",
        inherit.aes = FALSE
      )
  }

  .plt <- .plt + 
    facet_wrap(
      facets = vars({{facet_var}}),
      ncol = facet_ncol,
      scales = "free"
    ) + 
    labs(y = "proportion of all RL lineage cells") + 
    theme_light()
  
  return(.plt)
}
