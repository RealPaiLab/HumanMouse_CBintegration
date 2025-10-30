#' @import Seurat
#' @import qs

#' Subsets datasets to only RL lineage.
#' 
#' @param dataset_srat_path String indicating the path to dataset.
#' @param out_directory String indicating the path to write subsetted dataset.
#' @param author_name String indicating the authors name.
#' @param species String indicating species.
#' 



isolate_RL_lineage <- function(
    dataset_srat_path,
    out_directory,
    author_name,
    species
){
    # Defining RL lineage cells based on cell_type_mapping.csv
    rl_lineage = c("RL", "UBC", "UBC/GCP progenitor", "GCP", "GN", "oligodendrocyte/OPC", "microglia", "endothelial", "transitional cerebellar progenitor")  

    dataset_name <- paste(author_name, "RL", species, sep = "_")
    # Loading in dataset
    dataset_srat <- qread(dataset_srat_path)

    # Adding RL dataset metadata
    dataset_srat <- AddMetaData(object = dataset_srat, metadata = as.factor(dataset_name), col.name = "dataset_name")

    # Subsetting based on RL lineage
    rl_dataset_srat <- dataset_srat@meta.data[dataset_srat$common_cell_name %in% rl_lineage, ]
    rl_dataset <- subset(dataset_srat, cells = rownames(rl_dataset_srat))

    # Adding RL dataset metadata
    rl_dataset <- AddMetaData(object = rl_dataset, metadata = as.factor(dataset_name), col.name = "dataset_name")

    # Write subsetted dataset (qs)
    file_name_qs = paste(dataset_name, format(Sys.Date(), "%Y%m%d"), "stand.qs", sep = "_")
    qsave(rl_dataset, file.path(out_directory, file_name_qs))

    # Write subsetted dataset (rds)
    file_name_rds = paste(dataset_name, format(Sys.Date(), "%Y%m%d"), "stand.rds", sep = "_")
    saveRDS(rl_dataset, file.path(out_directory, file_name_rds))
}