# ==============================================================================
# Adds the Vladoiu annotations to the Seurat object
# ==============================================================================

# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# library(tidyverse)
# library(Seurat)

#' Adds Vladoiu annotations Seurat object metadata in the "cell_type" column
#' 
#' @param object Seurat object
#' 
add_annotations <- function(
  object
) {
  # set data and output directories
  root_dir <- "CBL_scRNAseq"
  data_dir <- file.path("/isilon", root_dir, "data/mouse/Vladoiu")
  
  # load data
  srat_annot <- readRDS(file.path(data_dir, "seurat_annotated.RDS"))
  annots <- srat_annot@active.ident %>% 
    data.frame(barcode = names(.), cell_type = ., row.names = NULL)
  
  # annotations are coded as E18.5 --> change to E18
  annots$barcode <- str_replace(annots$barcode, "E18.5", "E18")
  
  # match annotations using cell barcodes
  # merge on barcode
  metadata_with_annots <- left_join(
    x = object@meta.data,
    y = annots,
    by = "barcode"
  ) %>% 
    `rownames<-`(value = .$barcode)
  
  # replace metadata dataframe
  object@meta.data <- metadata_with_annots
  
  return(object)
}


#' Labels cell types using Vladoiu et al. 2019 Extended Data Figure 2.
#'
#' @param object Seurat object
#' @param metadata_col Column containing the unlabelled mouse cell types
#'
#' @return Seurat object with labelled cell types in the "mouse_cell_type"
#'   column.
#'
label_vladoiu_cells <- function(
  object, 
  metadata_col = "cell_type"
) {
  # add cell type labels to new column
  object@meta.data <- object@meta.data %>%
    mutate(
      mouse_cell_type = case_when(
        get(metadata_col) == "0" ~ "Excitatory cerebellar nuclei neurons",
        get(metadata_col) == "1" ~ "Embryonic and postnatal GCPs-1",
        get(metadata_col) == "2" ~ "Neural stem cells",
        get(metadata_col) == "3" ~ "Unipolar brush cell and GCP progenitor",
        get(metadata_col) == "4" ~ "Unipolar brush cells",
        get(metadata_col) == "5" ~ "GABA interneurons",
        get(metadata_col) == "6" ~ "Brainstem progenitors",
        get(metadata_col) == "7" ~ "Granule cells",
        get(metadata_col) == "8" ~ "VZ progenitors",
        get(metadata_col) == "9" ~ "Unipolar brush cell precursors",
        get(metadata_col) == "10" ~ "Differentiating Purkinje cells",
        get(metadata_col) == "11" ~ "Gliogenic progenitors-1",
        get(metadata_col) == "12" ~ "Upper rhombic lip progenitors",
        get(metadata_col) == "13" ~ "Mesenchymal stem cells-1",
        get(metadata_col) == "14" ~ "Purkinje cells",
        get(metadata_col) == "15" ~ "Postnatal GCPs-2",
        get(metadata_col) == "16" ~ "Post mitotic NTZ neurons",
        get(metadata_col) == "17" ~ "Roof plate-like stem cells",
        get(metadata_col) == "18" ~ "Proliferating VZ progenitors",
        get(metadata_col) == "19" ~ "Oligodendrocyte precursor cells",
        get(metadata_col) == "20" ~ "Gliogenic progenitors-2",
        get(metadata_col) == "21" ~ "Astrocyte/Bergmann glia precursors",
        get(metadata_col) == "22" ~ "Endothelial cells",
        get(metadata_col) == "23" ~ "Postnatal excitatory cerebellar nuclei neurons",
        get(metadata_col) == "24" ~ "GABA interneuron precursors",
        get(metadata_col) == "25" ~ "Pericytes",
        get(metadata_col) == "26" ~ "Early proliferating VZ progenitors",
        get(metadata_col) == "27" ~ "Mesenchymal stem cells-2",
        get(metadata_col) == "28" ~ "Microglia",
        get(metadata_col) == "29" ~ "Meninges",
        get(metadata_col) == "30" ~ "Red blood cells",
        TRUE ~ as.character(NA) # needs as.character, otherwise you get an error cuz data type isn't the same
      )
    )
  
  # change mouse NA cell types to "missing"
  if ("species" %in% colnames(object@meta.data)) {
    # if species are labelled
    object@meta.data <- object@meta.data %>%
      mutate(
        mouse_cell_type = case_when(
          species == "mouse" & is.na(get(metadata_col)) ~ "NA - missing",
          TRUE ~ as.character(mouse_cell_type)
        )
      )
  } else {
    # otherwise assume all cells are mouse cells
    object@meta.data <- object@meta.data %>%
      mutate(
        mouse_cell_type = case_when(is.na(get(metadata_col)) ~ "NA - missing",
                                    TRUE ~ as.character(mouse_cell_type)
        )
      )
  }
  
  # change cell type to factor and relevel by cell type lineage (e.g. UBCs,
  # granule cells, glial cells, etc.)
  object$mouse_cell_type <- factor(
    object$mouse_cell_type,
    levels = c(
      "Upper rhombic lip progenitors",
      "Unipolar brush cell and GCP progenitor",
      "Unipolar brush cell precursors",
      "Unipolar brush cells",
      "Embryonic and postnatal GCPs-1",
      "Postnatal GCPs-2",
      "Granule cells",
      "Post mitotic NTZ neurons",
      "Excitatory cerebellar nuclei neurons",
      "Postnatal excitatory cerebellar nuclei neurons",
      "Neural stem cells",
      "Mesenchymal stem cells-1",
      "Mesenchymal stem cells-2",
      "Roof plate-like stem cells",
      "Brainstem progenitors",
      "Early proliferating VZ progenitors",
      "Proliferating VZ progenitors",
      "VZ progenitors",
      "Differentiating Purkinje cells",
      "Purkinje cells",
      "GABA interneuron precursors",
      "GABA interneurons",
      "Gliogenic progenitors-1",
      "Gliogenic progenitors-2",
      "Astrocyte/Bergmann glia precursors",
      "Oligodendrocyte precursor cells",
      "Endothelial cells",
      "Meninges",
      "Microglia",
      "Pericytes",
      "Red blood cells",
      "NA - missing"
    )
  )
  
  return(object)
}
