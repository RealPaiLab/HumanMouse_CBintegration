# ==============================================================================
# functions to load data from Hendrikse et al. 2022
# ==============================================================================


#' Loads the mutated genes from Supplementary Table 5 of Hendrikse et al.
#' (Nature 2022). More details are in notes from 2023-02-16.
#'
#' @param mut_xlsx Path to file
#' @param rename_cols Whether or not to rename the columns
#'
#' @return Tibble of all mutated genes and p-values
#' 
load_mut_genes <- function(
  mut_xlsx = "/isilon/CBL_scRNAseq-archived/data/human/Hendrikse/supp_tab_5.xlsx",
  rename_cols = TRUE
) {
  mut_genes <- readxl::read_excel(
    path = mut_xlsx,
    sheet = 2,
    skip = 1
  )
  
  if (rename_cols) {
    mut_genes <- dplyr::rename(
      mut_genes,
      gene = Gene,
      freq_g4 = Frequency...3,
      pval_oncodrive_g4 = pval_OncoDriveFML...4,
      qval_oncodrive_g4 = qval_OncoDriveFML...5,
      pval_mutsig_g4 = pval_MutSig...6,
      qval_mutsig_g4 = qval_MutSig...7,
      freq_g3 = Frequency...8,
      pval_oncodrive_g3 = pval_OncoDriveFML...9,
      qval_oncodrive_g3 = qval_OncoDriveFML...10,
      pval_mutsig_g3 = pval_MutSig...11,
      qval_mutsig_g3 = qval_MutSig...12,
    )
  }
  
  return(mut_genes)
}
