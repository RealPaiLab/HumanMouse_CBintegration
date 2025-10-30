#' Convert a list of genes to its orthologous genes.
#' 
#' @param orig_genes Full list of genes which you want to convert from.
#' @param new_genes Full list of genes which you want to convert to.
#' @param gene_list List of genes that you want to convert.
#' 
get_orth_genes <- function(
  orig_genes,
  new_genes,
  gene_list
) {
  # convert old gene names to new gene names using the genes passed into the function
  orth_genes <- sapply(
    X = gene_list, 
    FUN = function(X) {
      if (X %in% orig_genes) {
        # if an orthologous gene exists, then set the name to the new gene
        X <- new_genes[orig_genes == X]
      } else {
        # if no orthologous gene was found, keep the old (original) gene name
        X
      }
    }, 
    USE.NAMES = FALSE
  )
  
  return(orth_genes)
}


#' Rename all genes in a Seurat object to their orthologs.
#' 
#' @param object Seurat object in which to rename the genes.
#' @param old_genes Vector of orthologs of the genes to rename from.
#' @param new_genes Vector of orthologs of the genes to rename to.
#' 
#' @return Seurat object with the renamed genes.
#' 
rename_genes <- function(
  object,
  old_genes,
  new_genes
) {
  # loop through assays
  for (assay in names(object@assays)) {
    
    # loop through counts, data, and scale.data slots
    for (slot in c("counts", "data", "scale.data")) {
      
      # get vector of old gene names
      old_names <- rownames(slot(object[[assay]], slot))
      
      # continue to next loop iteration if slot is empty
      if (length(old_names) == 0) {next}
      
      # get the new gene names
      new_names <- get_orth_genes(old_genes, new_genes, old_names)
      
      # set the new gene names
      rownames(slot(object[[assay]], slot)) <- unlist(new_names)
    }
    
    # set new names for var.features slot
    if (length(object[[assay]]@var.features) == 0) {
      next
    }
    old_names <- object[[assay]]@var.features
    object[[assay]]@var.features <- unlist(get_orth_genes(old_genes, new_genes, old_names))
  }
  
  return(object)
}


#' @import Seurat

#' Converts mouse genes to 1 to 1 orthologous human genes
#' 
#' @param dataset_srat Seurat object of dataset
#' @param ortho_gene_directory String that specifies location of orthologous gene


mouse_to_human_genes <- function(
  dataset_srat,
  ortho_gene_directory = "/u/llau/software/mb_scrnaseq/MB_scRNAseq/scrnaseq_Leo/utilities"
) {
  # Find file of orthologous genes 
  files <- list.files(ortho_gene_directory)
  pattern <- ".*hgnc_mgi_orth_genes.csv"
  matching_files <- grep(pattern, files, value = TRUE)

  if (length(matching_files) > 0) {
    message("Using existing orthologous genes file: ", matching_files[[1]])
    m2h_genes <- read.csv(file.path(ortho_gene_directory, matching_files[[1]]))
  } else {
    stop("No existing orthologous genes file")
  }

  # Remove all duplicated genes
  # Also remove Pisd (some wonky stuff is happening with that gene)
  m2h_genes <- m2h_genes[!(duplicated(m2h_genes$MGI.symbol) | duplicated(m2h_genes$MGI.symbol, fromLast = TRUE))
                        &!(duplicated(m2h_genes$HGNC.symbol) | duplicated(m2h_genes$HGNC.symbol, fromLast = TRUE))
                        & m2h_genes$MGI.symbol != "Pisd", ]

  # Change mouse gene names to human gene names
  # If human ortholog does not exist, keep mouse gene name
  dataset_srat_hgene <- rename_genes(dataset_srat, old_genes = m2h_genes$MGI.symbol, new_genes = m2h_genes$HGNC.symbol)
  return(dataset_srat_hgene)
}