
#' Plots UMAP for cells in each cluster, highlighting cells in red
#' @param obj (Seurat) input object
#' @param clName (char) metadata column on which to subset cells.
#' @return (ggarrange) multi-panel, one per cluster
DimPlot_HighlightSingleClusters <- function(obj, clName){

}

### FeaturePlot helper functions.
### From https://divingintogeneticsandgenomics.com/post/how-to-deal-with-overplotting-without-being-fooled/
##matrix_to_expression_df<- function(x, obj){
##        df<- x %>%
##                as.matrix() %>% 
##                as.data.frame() %>%
##                tibble::rownames_to_column(var= "gene") %>%
##                tidyr::pivot_longer(cols = -1, names_to = "cell", values_to = "expression") %>%
##                tidyr::pivot_wider(names_from = "gene", values_from = expression) %>%
##                left_join(obj@meta.data %>% 
##                                  tibble::rownames_to_column(var = "cell"))
##        return(df)
##}
##
##
##get_expression_data<- function(obj, assay = "RNA", slot = "data", 
##                               genes = NULL, cells = NULL){
##        if (is.null(genes) & !is.null(cells)){
##                df<- GetAssayData(obj, assay = assay, slot = slot)[, cells, drop = FALSE] %>%
##                        matrix_to_expression_df(obj = obj)
##        } else if (!is.null(genes) & is.null(cells)){
##                df <- GetAssayData(obj, assay = assay, slot = slot)[genes, , drop = FALSE] %>%
##                        matrix_to_expression_df(obj = obj)
##        } else if (is.null(genes & is.null(cells))){
##                df <- GetAssayData(obj, assay = assay, slot = slot)[, , drop = FALSE] %>%
##                        matrix_to_expression_df(obj = obj)
##        } else {
##                df<- GetAssayData(obj, assay = assay, slot = slot)[genes, cells, drop = FALSE] %>%
##                        matrix_to_expression_df(obj = obj)
##        }
##        return(df)
##}
