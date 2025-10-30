#' @import Seurat

#' Transfer missing labels
#' 
#' @param reference Seurat object of cells with annotations
#' @param query Seurat object of cells without annotations
#' @param predict_column Name of column containing annotation "cell type name"


label_transfer <- function(
  reference,
  query,
  predict_column
) {
    anchors <- FindTransferAnchors(reference = reference, query = query)
    predictions <- TransferData(anchorset = anchors, refdata = reference@meta.data[[predict_column]])
    query@meta.data[[predict_column]] <- predictions$predicted.id
    dataset_srat<- merge(x = reference, y = query, add.cell.ids = c("Reference", "Query"))


    reference@meta.data$reference_cell_type <- reference@meta.data[[predict_column]] 
    query@meta.data$query_cell_type <- query@meta.data[[predict_column]]
    query <- AddMetaData(object = query, metadata = predictions) 
    dataset_srat_csv<- merge(x = reference, y = query, add.cell.ids = c("Reference", "Query"))
    metadata <- as.data.frame(dataset_srat_csv@meta.data)

    label_list <- list(dataset_srat, metadata)
    return(label_list)
}

