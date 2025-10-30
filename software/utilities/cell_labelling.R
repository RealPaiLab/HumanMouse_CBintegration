# ==============================================================================
# Functions for labelling cells and cell types.
# ==============================================================================


#' Combine Aldinger and Vladoiu cell types to a single column with common names.
#' 
#' @param srat Seurat object containing labelled human cells and mouse cells.
#' @param col_name New column name for the pooled cell types.
#' 
#' @return Seurat object with the pooled cell types in the metadata.
#' 
pool_cell_types <- function(
  srat,
  col_name = "common_cell_type"
) {
  if ("new_cell_type" %in% colnames(srat@meta.data)) {
    # `new_cell_type` contains Liam's annotations; use these annotations if they're available
    case_aldinger <- case_aldinger_liam
    message("Using the `new_cell_type` metadata column to generate a common cell type")
  } else if ("figure_clusters" %in% colnames(srat@meta.data)) {
    # `figure_clusters` contains the original Aldinger annotations
    case_aldinger <- case_aldinger_orig
    message("Using the `figure_clusters` metadata column to generate a common cell type")
  }
  
  srat@meta.data <- mutate(
    srat@meta.data,
    aldinger = case_aldinger(srat),
    vladoiu = case_vladoiu(srat),
    {{col_name}} := coalesce(aldinger, vladoiu) %>% replace_na("other/missing")
  ) %>% 
    select(!c(aldinger, vladoiu))
  
  return(srat)
}


#' For internal use only with the `pool_cell_types` function. Converts Liam's
#' annotations of the Aldinger data into a common annotation.
#'
#' @param srat Seurat object containing labelled human cells and mouse cells.
#'
#' @return A `case_when` pattern for use in the `mutate` function.
#' 
case_aldinger_liam <- function(
  srat
) {
  new_cell_type <- srat$new_cell_type
  case_when(
    new_cell_type %in% c("RL-VZ", "RL-SVZ") ~ "RL",
    new_cell_type %in% c("Early UBCs", "Late UBCs") ~ "UBCs",
    new_cell_type %in% c("GCP") ~ "GCPs",
    new_cell_type %in% c("Early GN", "GN") ~ "GNs",
    TRUE ~ NA_character_
  )
}


#' For internal use only with the `pool_cell_types` function. Converts the
#' original Aldinger annotations into a common annotation.
#'
#' @param srat Seurat object containing labelled human cells and mouse cells.
#'
#' @return A `case_when` pattern for use in the `mutate` function.
#' 
case_aldinger_orig <- function(
  srat
) {
  cell_type <- srat$figure_clusters
  case_when(
    cell_type %in% c("01-PC") ~ "PCs",
    cell_type %in% c("02-RL") ~ "RL",
    cell_type %in% c("03-GCP") ~ "GCPs",
    cell_type %in% c("04-GN") ~ "GNs",
    cell_type %in% c("05-eCN/UBC") ~ "UBCs",
    cell_type %in% c("06-iCN") ~ "cerebellar nuclei",
    cell_type %in% c("07-PIP", "18-MLI") ~ "interneurons",
    cell_type %in% c("08-BG", "10-Glia") ~ "glia",
    cell_type %in% c("09-Ast", "19-Ast/Ependymal") ~ "astrocytes",
    cell_type %in% c("11-OPC", "12-Committed OPC") ~ "OPCs",
    cell_type %in% c("13-Endothelial") ~ "endothelial",
    cell_type %in% c("14-Microglia") ~ "microglia",
    cell_type %in% c("15-Meninges") ~ "meninges",
    cell_type %in% c("16-Pericytes") ~ "pericytes",
    cell_type %in% c("17-Brainstem") ~ "brainstem",
    cell_type %in% c("20-Choroid", "21-BS Choroid/Ependymal") ~ "choroid",
    TRUE ~ NA_character_
  )
}


#' For internal use only with the `pool_cell_types` function. Converts the
#' original Vladoiu annotations into a common annotation.
#'
#' @param srat Seurat object containing labelled human cells and mouse cells.
#'
#' @return A `case_when` pattern for use in the `mutate` function.
#' 
case_vladoiu <- function(
  srat
) {
  mouse_cell_type <- srat$mouse_cell_type
  case_when(
    mouse_cell_type %in% c(
      "Unipolar brush cell and GCP progenitor",
      "Unipolar brush cell precursors"
    ) ~ "RL",
    mouse_cell_type %in% c("Unipolar brush cells") ~ "UBCs",
    mouse_cell_type %in% c(
      "Embryonic and postnatal GCPs-1",
      "Postnatal GCPs-2"
    ) ~ "GCPs",
    mouse_cell_type %in% c("Granule cells") ~ "GNs",
    mouse_cell_type %in% c(
      "Excitatory cerebellar nuclei neurons",
      "Postnatal excitatory cerebellar nuclei neurons"
    ) ~ "cerebellar nuclei",
    # mouse_cell_type %in% c(
    #   "Upper rhombic lip progenitors",
    #   "Neural stem cells",
    #   "Mesenchymal stem cells-1",
    #   "Mesenchymal stem cells-2",
    #   "Roof plate-like stem cells",
    #   "Brainstem progenitors"
    # ) ~ "stem cells",
    mouse_cell_type %in% c(
      "Early proliferating VZ progenitors",
      "Proliferating VZ progenitors",
      "VZ progenitors"
    ) ~ "VZ progenitors",
    mouse_cell_type %in% c(
      "Differentiating Purkinje cells",
      "Purkinje cells"
    ) ~ "PCs",
    mouse_cell_type %in% c(
      "GABA interneuron precursors",
      "GABA interneurons"
    ) ~ "interneurons",
    mouse_cell_type %in% c(
      "Gliogenic progenitors-1",
      "Gliogenic progenitors-2"
    ) ~ "glia",
    mouse_cell_type %in% c("Astrocyte/Bergmann glia precursors") ~ "astrocytes",
    mouse_cell_type %in% c("Oligodendrocyte precursor cells") ~ "OPCs",
    mouse_cell_type %in% c("Microglia") ~ "microglia",
    mouse_cell_type %in% c("Endothelial cells") ~ "endothelial",
    mouse_cell_type %in% c("Meninges") ~ "meninges",
    mouse_cell_type %in% c("Pericytes") ~ "pericytes",
    TRUE ~ NA_character_
  )
}


#' Label cells with TRUE/FALSE based on clustering or other cell metadata.
#'
#' For example, the function can be used to label human-specific UBC cells based
#' on the `seurat_clusters` metadata column by specifying which clusters are
#' human-specific and which are not.
#'
#' @param srat Seurat object.
#' @param eval_value Value of `eval_name` to label as TRUE.
#' @param eval_name Column to evaluate for returning TRUE/FALSE.
#' @param col_name New column name for cell metadata.
#' 
#' @return Seurat object with the new cell metadata column.
#' 
add_cell_label <- function(
  srat,
  eval_value,
  eval_name = "seurat_clusters",
  col_name = "human_specific"
) {
  #TODO
}


#' Label cell type of integrated clusters.
#'
#' @param srat Seurat object
#' @param col_name New column name for the cell type labels
#' @param integ_method One of `c("cca", "harmony")`
#' @param datasets One of `c("rl", "full")` where "rl" is for integration of RL
#'   only and "full" is for full integration of Aldinger and Vladoiu datasets.
#'
#' @return
#' @export
#'
#' @examples
label_integ_clusters <- function(
  srat,
  col_name = "integ_clusters",
  integ_method = "harmony",
  datasets = "rl"
) {
  # function to use for labelling clusters
  case_cluster <- switch(
    integ_method,
    "cca" = switch(datasets, "rl" = case_cca, "full" = case_cca_full),
    "harmony" = case_harmony,
    stop('`integ_method` must be one of c("cca", "harmony")')
  )
  
  # label clusters with cell type
  srat@meta.data <- mutate(
    srat@meta.data,
    {{col_name}} := case_cluster(srat) %>% as.factor
  )
  
  # relevel cluster names in numerical order
  srat@meta.data[[col_name]] <- fct_relevel(
    srat@meta.data[[col_name]],
    ~ str_sort(levels(srat@meta.data[[col_name]]), numeric = TRUE)
  )
  
  return(srat)
}


#' For internal use only with the `label_integ_clusters` function. Adds cell
#' types to the Seurat clusters (based on the published annotations).
#'
#' This function should only be used on "<results>/vladoiu_liam_RL.rds".
#'
#' @param srat Seurat object.
#'
#' @return A `case_when` pattern for use in the `mutate` function.
#' 
case_cca <- function(
  srat
) {
  cluster_num <- srat$seurat_clusters
  
  # see:
  # <results>/20220911/num_cells_per_cluster.pdf
  # <results>/20220911/human_cells_per_cluster.pdf
  # <results>/20220916/mouse_cells_per_cluster.pdf
  case_when(
    cluster_num == "0" ~ "0-Mouse GCP",
    cluster_num == "1" ~ "1-Mouse GCP",
    cluster_num == "2" ~ "2-Mouse GN",
    cluster_num == "3" ~ "3-Mouse GN",
    cluster_num == "4" ~ "4-Mouse GCP",
    cluster_num == "5" ~ "5-Human GN and Mouse UBC",
    cluster_num == "6" ~ "6-Human GN",
    cluster_num == "7" ~ "7-Homol UBC",
    cluster_num == "8" ~ "8-Human GCP/RL_SVZ",
    cluster_num == "9" ~ "9-Mouse UBC",
    cluster_num == "10" ~ "10-Human GN/Mouse UBC",
    cluster_num == "11" ~ "11-Mouse UBC Precursor",
    cluster_num == "12" ~ "12-Human RL_VZ and Mouse GCP/progenitor",
    cluster_num == "13" ~ "13-Mouse GCP",
    cluster_num == "14" ~ "14-Mouse GCP/progenitor",
    cluster_num == "15" ~ "15-Unknown",
    cluster_num == "16" ~ "16-Mouse GCP/progenitor",
    cluster_num == "17" ~ "17-Unknown",
    cluster_num == "18" ~ "18-Mouse GN",
    cluster_num == "19" ~ "19-NonHomol UBC",
    cluster_num == "20" ~ "20-NonHomol UBC",
    cluster_num == "21" ~ "21-Human GN",
    cluster_num == "22" ~ "22-Unknown"
  )
}


#' For internal use only with the `label_integ_clusters` function. Adds cell
#' types to the Seurat clusters (based on the published annotations).
#' 
#' This function should only be used on "<results>/20230126/without_future/aldinger_vladoiu_cca.rds".
#'
#' @param srat Seurat object.
#'
#' @return A `case_when` pattern for use in the `mutate` function.
#'
# case_cca_full <- function(
#   srat
# ) {
#   cluster_num <- srat$seurat_clusters
#   
#   # see:
#   # <results>/20230129/num_cells_per_cluster.pdf
#   # <results>/20230129/human_cells_per_cluster.pdf
#   # <results>/20230129/mouse_cells_per_cluster.pdf
#   case_when(
#     cluster_num == "0" ~ "PCs",
#     cluster_num == "1" ~ "Cerebellar nuclei/MLI",
#     cluster_num == "2" ~ "Glia",
#     cluster_num == "3" ~ "UBCs/GNs",
#     cluster_num == "4" ~ "GCPs",
#     cluster_num == "5" ~ "Astrocytes",
#     cluster_num == "6" ~ "PCs",
#     cluster_num == "7" ~ "GCPs",
#     cluster_num == "8" ~ "Interneurons",
#     cluster_num == "9" ~ "PCs",
#     cluster_num == "10" ~ "UBCs",
#     cluster_num == "11" ~ "PCs",
#     cluster_num == "12" ~ "GNs",
#     cluster_num == "13" ~ "Cerebellar nuclei",
#     cluster_num == "14" ~ "NSCs",
#     cluster_num == "15" ~ "Interneurons",
#     cluster_num == "16" ~ "UBCs/GNs",
#     cluster_num == "17" ~ "BS/choroid",
#     cluster_num == "18" ~ "GNs",
#     cluster_num == "19" ~ "Progenitors",
#     cluster_num == "20" ~ "GCPs",
#     cluster_num == "21" ~ "Glia",
#     cluster_num == "22" ~ "Choroid",
#     cluster_num == "23" ~ "PCs",
#     cluster_num == "24" ~ "Interneurons",
#     cluster_num == "25" ~ "OPCs",
#     cluster_num == "26" ~ "",
#     cluster_num == "27" ~ "",
#     cluster_num == "28" ~ "",
#     cluster_num == "29" ~ "",
#     cluster_num == "30" ~ "",
#     cluster_num == "31" ~ "",
#     cluster_num == "32" ~ "",
#     cluster_num == "33" ~ "",
#     cluster_num == "34" ~ "",
#     cluster_num == "35" ~ "",
#     cluster_num == "36" ~ "",
#     cluster_num == "37" ~ "OPCs",
#     cluster_num == "38" ~ "",
#     cluster_num == "39" ~ "",
#     cluster_num == "40" ~ "",
#   )
# }


#' For internal use only with the `label_integ_clusters` function. Adds cell
#' types to the Seurat clusters (based on the published annotations).
#'
#' This function should only be used on "<results>/20221003/vladoiu_liam_RL_harmony.rds".
#'
#' @param srat Seurat object.
#'
#' @return A `case_when` pattern for use in the `mutate` function.
#' 
case_harmony <- function(
  srat
) {
  cluster_num <- srat$seurat_clusters
  case_when(
    cluster_num == "0" ~ "0-GN",
    cluster_num == "1" ~ "1-GCP",
    cluster_num == "2" ~ "2-Early GN",
    cluster_num == "3" ~ "3-GCP",
    cluster_num == "4" ~ "4-GCP",
    cluster_num == "5" ~ "5-Homol UBC",
    cluster_num == "6" ~ "6-GCP",
    cluster_num == "7" ~ "7-Mouse UBC",
    cluster_num == "8" ~ "8-GN",
    cluster_num == "9" ~ "9-UBC Precursor",
    cluster_num == "10" ~ "10-GN",
    cluster_num == "11" ~ "11-RL_SVZ",
    cluster_num == "12" ~ "12-GCP",
    cluster_num == "13" ~ "13-GCP",
    cluster_num == "14" ~ "14-RL_VZ and GCP",
    cluster_num == "15" ~ "15-Early GN",
    cluster_num == "16" ~ "16-NonHomol UBC",
    cluster_num == "17" ~ "17-GN",
    cluster_num == "18" ~ "18-GCP",
    cluster_num == "19" ~ "19-GN",
    cluster_num == "20" ~ "20-GCP",
    cluster_num == "21" ~ "21-Mouse UBC",
    cluster_num == "22" ~ "22-Unknown",
    cluster_num == "23" ~ "23-GCP",
    cluster_num == "24" ~ "24-Unknown"
  )
}


#' Label cells from the RL lineage integration.
#' 
#' This function should only be used on `<llau integrated
#' results>/20240524/25_pc_without_luo/25_pc_rl.qs`.
#'
#' @param srat Seurat object in <llau integrated
#'   results>/20240524/25_pc_without_luo/25_pc_rl.qs
#'
#' @return Seurat object with the UBC subclusters (`ubc_subclusters` and
#'   `final_clusters`) and cell type annotations (`cell_type_annot`) in the
#'   metadata.
#'
label_rl_lineage_integration <- function(
  srat
) {
  # load Seurat metadata containing UBC subclusters
  ubc_subclusters <- read.csv(
    "/.mounts/labs/pailab/private/llau/results/integrated/20240527/ubc_subset_metadata.csv",
    row.names = 1 # row names in the first column
  ) %>%
    dplyr::rename(ubc_subclusters = snn_res.0.3) %>%
    dplyr::select(ubc_subclusters)

  # add UBC subclusters to the original Seurat metadata
  md <- merge(
    x = srat[[]],
    y = ubc_subclusters,
    by = "row.names",
    all.x = TRUE
  ) %>%
    mutate(
      # combine UBC subclusters with other clusters
      final_clusters = case_when(
        !is.na(ubc_subclusters) ~ paste0("UBC_", ubc_subclusters),
        .default = snn_res.0.4
      ),
      # set factor levels for the new column
      final_clusters = fct_relevel(
        final_clusters,
        str_sort(names(table(final_clusters)), numeric = TRUE)
      )
    ) %>%
    column_to_rownames(var = "Row.names")

  # load cell type annotations
  cluster_annot <- read_csv("/.mounts/labs/pailab/private/llau/results/integrated/20240618/clusters_cell_types.csv") %>%
    mutate(
      final_clusters = str_remove(Cluster, "^Cluster_"),
      .before = 1
    ) %>%
    dplyr::rename(cell_type_annot = broad_w_ubc_subtypes) %>%
    dplyr::select(final_clusters, cell_type_annot)

  # add cell type annotations to the metadata
  md <- merge(
    x = rownames_to_column(md, var = "cell_id"),
    y = cluster_annot,
    by = "final_clusters",
    all.x = TRUE
  ) %>%
    relocate(final_clusters, .after = ubc_subclusters) %>%
    column_to_rownames(var = "cell_id")

  # sort metadata rows back to original order and add back to the Seurat object
  md <- md[rownames(srat[[]]), ]
  srat@meta.data <- md

  return(srat)
}
