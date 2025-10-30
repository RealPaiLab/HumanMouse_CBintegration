#' @importFrom biomaRt useMart getLDS

# ==============================================================================
# see https://github.com/satijalab/seurat/issues/2493 and
# https://www.r-bloggers.com/2016/10/converting-mouse-to-human-gene-names-with-biomart-package/
# ==============================================================================

# load human biomart database
load_human_bm <- function() {
  host <- "https://dec2021.archive.ensembl.org"
  message(sprintf("Using Ensembl host at %s for human BioMart", host))
  return(useEnsembl(biomart = "ensembl", 
                    dataset = "hsapiens_gene_ensembl", 
                    host = host))
}

# load mouse biomart database
load_mouse_bm <- function() {
  host <- "https://dec2021.archive.ensembl.org"
  message(sprintf("Using Ensembl host at %s for mouse BioMart", host))
  return(useEnsembl(biomart = "ensembl", 
                    dataset = "mmusculus_gene_ensembl", 
                    host = host))
}

# Convert human genes to mouse genes
convert_h2m <- function(
  human_bm,
  mouse_bm,
  human_genes = NULL
) {
  # if genes aren't provided, then don't filter biomart request;
  # otherwise, only return the filtered genes
  if (is.null(human_genes)) {
    filters = ""
    values = ""
  } else {
    filters = "hgnc_symbol"
    values = human_genes
  }
  
  # get homologous genes
  mouse_genes <- getLDS(
    attributes = c("hgnc_symbol"),
    filters = filters,
    values = values,
    mart = human_bm,
    attributesL = c("mgi_symbol"),
    martL = mouse_bm,
    uniqueRows = TRUE
  )
  
  return(mouse_genes)
}

# Convert mouse genes to human genes
convert_m2h <- function(
  mouse_bm,
  human_bm,
  mouse_genes = NULL
) {
  # if genes aren't provided, then don't filter biomart request;
  # otherwise, only return the filtered genes
  if (is.null(mouse_genes)) {
    filters = ""
    values = ""
  } else {
    filters = "hgnc_symbol"
    values =  mouse_genes
  }
  
  # get homologous genes
  human_genes <- getLDS(
    attributes = c("mgi_symbol"),
    filters = filters,
    values = values,
    mart = mouse_bm,
    attributesL = c("hgnc_symbol"),
    martL = human_bm,
    uniqueRows = TRUE
  )
  
  return(human_genes)
}


#' Load orthologous genes from file and remove duplicates.
#'
#' @param file Path to CSV of orthologous genes.
#' @param filter_duplicates If `TRUE`, only genes that match one-to-one will be
#'   returned.
#' @param cnames Column names to check for non-one-to-one orthologs. Default is
#'   `MGI.symbol` and `HGNC.symbol`.
#'
#' @return A data frame with the orthologous genes.
#' 
load_orth_genes <- function(
  file,
  filter_duplicates = TRUE,
  cnames = c("MGI.symbol", "HGNC.symbol")
) {
  orth_gene_list <- readr::read_csv(file = file)
  nrow_start <- nrow(orth_gene_list)
  
  if (filter_duplicates) {
    # remove all duplicated genes
    # also remove Psid (some wonky stuff is happening with that gene)
    orth_gene_list <- dplyr::filter(
      orth_gene_list,
      if_all(all_of(cnames), ~ duplicated(.) == FALSE) & 
        if_any(all_of(cnames), ~ . != "Pisd", )
    )
  }
  nrow_end <- nrow(orth_gene_list)
  
  message(sprintf(
    "Filtered out %s genes with multiple orthologs",
    nrow_start - nrow_end
  ))
  
  return(orth_gene_list)
}


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


#' [DEPRECATED: use `get_orth_genes` instead.] Get orthologous gene names from
#' their original genes.
#' 
#' For example, to convert mouse genes to their orthologous human genes, you
#' need a table of matching orthologous genes, with the mouse gene column passed
#' to `old_genes` and the human gene column passed to `new_genes`. The list of
#' mouse genes to be converted is passed to `old_names`.
#' 
#' If the gene ortholog is missing from orthology table, the old gene name will
#' be kept.
#' 
#' @param old_genes Vector of orthologs of the genes to convert from.
#' @param new_genes Vector of orthologs of the genes to convert to.
#' @param old_names Vector of genes to convert from.
#' 
#' @return Vector of orthologous genes in the same order as the input argument.
#' 
get_new_names <- function(
  old_genes,
  new_genes,
  old_names
) {
  new_names <- get_orth_genes(old_genes, new_genes, old_names)
  return(new_names)
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
