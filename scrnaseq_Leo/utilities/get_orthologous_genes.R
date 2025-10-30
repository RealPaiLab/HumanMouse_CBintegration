#' @import biomaRt  
#' @import useMart
#' @import getLDS

#' Loads human biomart database


load_human_bm <- function() {
    return(useEnsembl(
        biomart = "ensembl", 
        dataset = "hsapiens_gene_ensembl", 
        host = "https://dec2021.archive.ensembl.org"))
}



#' Loads mouse biomart database


load_mouse_bm <- function() {
    return(useEnsembl(
        biomart = "ensembl", 
        dataset = "mmusculus_gene_ensembl", 
        host = "https://dec2021.archive.ensembl.org"))
}



#' Converts human genes to mouse genes


convert_h2m <- function(
    human_genes, 
    human_bm, 
    mouse_bm
) {
    # get homologous genes
    mouse_genes <- getLDS(attributes = c("hgnc_symbol"), 
                            filters = "hgnc_symbol", 
                            values = human_genes, 
                            mart = human_bm, 
                            attributesL = c("mgi_symbol"), 
                            martL = mouse_bm, 
                            uniqueRows = TRUE)
    
    return(mouse_genes)
}


#' @import Seurat
#' @import clustree

#' Gets orthologous mice genes
#' 
#' @param dataset_srat Seurat object containing mouse dataset
#' @param ortho_gene_directory String that specifies location to print orthologous genes csv file


get_orthologous_genes <- function(
    dataset_srat,
    ortho_gene_directory
) {
    #load biomart/ensembl databases
    message("Loading biomart/ensembl databases")
    human_bm <- load_human_bm()
    mouse_bm <- load_mouse_bm()

    # get mouse genes and convert to human genes
    message("Converting mouse to human genes")
    mouse_genes <- rownames(dataset_srat)
    m2h_genes <- convert_m2h(mouse_genes, mouse_bm, human_bm)
    
    # save the orthologous genes to a CSV file
    write.csv(m2h_genes, file = file.path(ortho_gene_directory, paste0(format(Sys.Date(), "%Y%m%d"),"_hgnc_mgi_orth_genes.csv")))
}


    