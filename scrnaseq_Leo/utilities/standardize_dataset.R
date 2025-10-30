#' @import Seurat
#' @import tidyverse

#' QC datasets
#' Filtering based on percent_mt and nFeatures
#' 
#' @param dataset_srat Seurat object of dataset
#' @param mito_pattern String that specifies mitochondrial gene identifier

# Setting working directory
setwd("/u/llau/software/mb_scrnaseq/MB_scRNAseq")

qc_dataset <- function(
  dataset_srat,
  mito_pattern = "^MT-"
) {
    # Search if percent_mt already exists. Can exist in many formats 
    matching_percent_mt <- grep("^percent[._]m.*", colnames(dataset_srat@meta.data), value = TRUE)

    # Determine mitochondrial percentage
    if(length(matching_percent_mt) == 1) {
      # Standardizing percent_mt name
      dataset_srat[["percent_mt"]] <- dataset_srat[[matching_percent_mt]]
      dataset_srat[[matching_percent_mt]] <- NULL
    } else {
      # Generates percent_mt column based on metadata
      dataset_srat[["percent_mt"]] <- PercentageFeatureSet(dataset_srat, pattern = mito_pattern)
    }

    # Defining QC metrics
    mt_limit <- median(dataset_srat$percent_mt) + (4 * sd(dataset_srat$percent_mt))
    feat_limit <- median(dataset_srat$nFeature_RNA) + (4 * sd(dataset_srat$nFeature_RNA))

    # Filter out low-quality cells (< 200 genes) and cell multiplets ( > 4 S.D) and mitochondiral counts (> 4 S.D)
    # 4 S.D metric used by Vladoiu 2019 paper
    dataset_srat <- subset(dataset_srat, subset = nFeature_RNA >= 200 & nFeature_RNA <= feat_limit & percent_mt <= mt_limit)
    
    return(dataset_srat)
}




#' @import Seurat
#' @import tidyverse
#' @import qs

#' Standardize datasets
#' Adding identifier metadata -> SCTransform and dimensionality reduction -> Plotting TSNE and UMAP -> Exporting standardized dataset
#' 
#' @param dataset_srat Seurat object of dataset
#' @param dataset_author String that specifies dataset author
#' @param dataset_name String that specifies dataset name
#' @param species String that specifies dataset species
#' @param regress_var String that specifies dataset column to regress
#' @param cell_type_column String that specifies dataset column containing cell type identifier
#' @param do_qc Boolean indicating whether or not to QC dataset
#' @param out_directory String that specifies output directory


standardize_dataset <- function(
  dataset_srat,
  dataset_author,
  dataset_name,
  species,
  regress_var = NULL,
  cell_type_column,
  do_qc = FALSE,
  out_directory
) {
    #QC dataset
    #if(do_qc){
      #dataset_srat = qc_dataset(dataset_srat = dataset_srat, mito_pattern = mito_pattern)
    #}

    # Adding species and dataset metadata
    dataset_srat <- AddMetaData(object = dataset_srat, metadata = as.factor(species), col.name = "species") %>%
      AddMetaData(metadata = as.factor(dataset_name), col.name = "dataset_name")
    
    # Renaming cell type column to standardize
    if(cell_type_column != "cell_type"){
      dataset_srat$cell_type <- dataset_srat[[cell_type_column]]
      dataset_srat[[cell_type_column]] <- NULL
    }
    
    # Adding common cell type metadata
    common_cell_types <- read.csv("scrnaseq_Leo/utilities/cell_type_mapping.csv")
    original_cell_names <- dataset_srat$cell_type
    matching_indices <- match(original_cell_names, common_cell_types$original)
    dataset_srat$common_cell_name <- common_cell_types$common[matching_indices]

    # Making sure this column is a factor
    dataset_srat$common_cell_name <- factor(dataset_srat$common_cell_name)

    # Adding author name to cell type
    dataset_srat$cell_type <- paste(dataset_srat$cell_type, dataset_author, sep = "_")

    # Removing previously generated clusters
    pattern <- paste0(".*_snn_res.")
    matching_column <- grep(pattern, colnames(dataset_srat@meta.data), value = TRUE)
    for(column in matching_column){
      dataset_srat[[column]] <- NULL
    }

    # Converting mouse to human genes
    source("scrnaseq_Leo/utilities/mouse_to_human_genes.R")
    if( species == "mouse"){
      message("Dataset ", dataset_name, " is mice, converting genes")
      dataset_srat <- mouse_to_human_genes(dataset_srat = dataset_srat)
    } else {
      message("Dataset ", dataset_name, " is human, will skip this step")
    }

    # Perform SCTransform and dimensionality reduction suite (PCA)
    dataset_srat <- SCTransform(dataset_srat, vars.to.regress = regress_var, variable.features.n = 5000, return.only.var.genes = FALSE) %>%
      RunPCA()

    # Save QC'd dataset as a .rds file
    file_name = paste(dataset_name, species, format(Sys.Date(), "%Y%m%d"), "stand.rds", sep = "_")
    saveRDS(object = dataset_srat, file = file.path(out_directory, file_name))

    # Save QC'd dataset as a .qs file
    file_name_qs = paste(dataset_name, species, format(Sys.Date(), "%Y%m%d"), "stand.qs", sep = "_")
    qsave(dataset_srat, file.path(out_directory, file_name_qs))

    message("Dataset successfully standardized. The standardized dataset is called: ", file_name)
}

