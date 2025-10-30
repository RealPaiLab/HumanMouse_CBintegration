#' @import Seurat
#' @import qs
#' @import tidyverse

#' Subsets datasets to only RL lineage.
#' 
#' @param integ_fc Seurat object.
#' 

fc_subset <- function(
    integ_fc
){
  integ_fc@meta.data <- integ_fc@meta.data %>% 
    mutate(
      # timepoints from human datasets
      human_age = case_when(
        dataset_name == "Aldinger_full_cerebellum_human" ~ age,
        dataset_name == "Luo_full_cerebellum_human" ~ paste(
          str_remove(string = batch, pattern = "PCW"),
          "PCW"
        ),
        dataset_name == "Sepp_full_cerebellum_human" ~ paste(
          str_remove(string = Stage, pattern = " wpc"),
          "PCW"
        )
      ),
      # timepoints from mouse datasets
      mouse_age = case_when(
        dataset_name == "Vladoiu_full_cerebellum_mouse" ~ str_remove(
          string = orig.ident,
          pattern = "Vladoiu-"
        ),
        dataset_name == "Sepp_full_cerebellum_mouse" ~ Stage
      )
    )
  integ_fc$combined_age <- ifelse(integ_fc$species == "human", integ_fc$human_age, integ_fc$mouse_age)

  # convert ages to factor
  integ_fc$human_age <- fct(
    x = integ_fc$human_age,
    levels = str_sort(unique(integ_fc$human_age), numeric = TRUE, na_last = NA)
  )
  integ_fc$mouse_age <- fct(
    x = integ_fc$mouse_age,
    levels = str_sort(unique(integ_fc$mouse_age), numeric = TRUE, na_last = NA)
  )
  integ_fc$combined_age <- fct(
    x = integ_fc$combined_age,
    levels = str_sort(unique(integ_fc$combined_age), numeric = TRUE, na_last = NA)
  )

  # Subsetting out cells younger than 11 PCW
  Idents(integ_fc) <- "combined_age"
  integ_fc_age_subset <- subset(x = integ_fc, idents = c("7 PCW", "8 PCW", "9 PCW", "10 PCW"), invert = TRUE)

  # Subset for RL cells 
  Idents(integ_fc_age_subset) <- "common_cell_name"
  rl_cells <- c("RL", "UBC", "UBC/GCP progenitor", "GCP", "GN", "oligodendrocyte/OPC", "microglia", "endothelial")  
  integ_rl <- subset(x = integ_fc_age_subset, idents = rl_cells)

  return(integ_rl)
}