#' @import Seurat

#' Prepare dataset for integration by appending datsets in a list, applying SCTransform() and RunPCA() to each
#' 
#' @param dataset_list List of Seurat object of datasets
#' @param run_sctransform Boolean indicating whether to run SCTransform before integration


integration_prep <- function(
  dataset_list,
  run_sctransform = FALSE
) {
  # re-run sctransform
  if(run_sctransform){
    dataset_list <- lapply(X = dataset_list, FUN = SCTransform,
                    variable.features.n = 5000,
                    return.only.var.genes = FALSE)
  }
  
  sex_genes <- c(
  "XIST", "DDX3Y", "KDM5D", "RPS4Y1", "USP9Y", "UTY", # Vawter et al. 2003
  "Xist", "Tsix", "Eif2s3y", "Ddx3y", "Uty", "Kdm5d", # Hamed et al. 2022
  "TSIX", "EIF2S3B", # human orthologs of Hamed mouse genes
  "Usp9y" # mouse ortholog of Vawter human genes
  )
  
  dataset_list <- lapply(
    X = dataset_list,
    function(obj) {
    subset(obj, features = setdiff(rownames(obj),sex_genes))
    }
  )

  # prep for integration
  features <- SelectIntegrationFeatures(dataset_list, nfeatures = 5000)

  dataset_list <- lapply(
    X = dataset_list,
    FUN = RunPCA,
    features = features,
    npcs = 100
  )

  dataset_list <- PrepSCTIntegration(dataset_list, anchor.features = features)

  integration_vars <- list(dataset_list = dataset_list, integ_feat = features)
  return(integration_vars)
}