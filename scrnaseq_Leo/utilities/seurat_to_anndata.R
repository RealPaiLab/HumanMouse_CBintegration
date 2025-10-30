#' @import Seurat
#' @import SeuratDisk

#' Convert Seurat object into Anndata object
#' 
#' @param srat_object Seurat object containing dataset
#' @param assay String or vector of name(s) of assay(s) to keep
#' @param dimreducs String or vector of name(s) of dimensionality reduction(s) to keep
#' @param out_name String that specifies output file name
#' @param out_directory String that specifies output directory


seurat_to_anndata <- function(
  srat_object,
  assay,
  dimreducs,
  out_name,
  out_directory
) {
    srat_filtered <- DietSeurat(srat_object, counts = T, scale.data = F, features = NULL,
           assays = assay, dimreducs = dimreducs, graphs = NULL)

    for (col in colnames(srat_filtered@meta.data)[sapply(srat_filtered@meta.data, class) == "factor"]) {
      srat_filtered[[col]] <- as.character(srat_filtered@meta.data[[col]])
    }

    file_name <- paste0(out_name, ".loom")
    SaveLoom(srat_filtered, filename = paste(output_path, file_name, sep = "/"))
}
