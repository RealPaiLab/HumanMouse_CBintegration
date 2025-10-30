#' @import qs
#' @import future

#' Coverting .rds files to .qs files for faster read/write.
#' 
#' @param rds_dataset_path String indicating path to .rds dataset


rds_to_qs <- function(
  rds_dataset_path
) {
    srat_dataset <- readRDS(rds_dataset_path)
    # Replace ".rds" with ".qs"
    qs_dataset_path <- gsub("\\.rds$", ".qs", rds_dataset_path)
    # Save dataset as .qs file
    qsave(srat_dataset, qs_dataset_path)
}

library(qs)
library(future)

setwd("/u/llau/software/mb_scrnaseq/MB_scRNAseq")
datasets <- read.csv("scrnaseq_Leo/integrations/dataset_list.csv")

plan(multicore)

futures <- list()

for(i in seq_along(datasets)){
  futures[[i]] <- future({
      rds_to_qs(datasets[i, "dataset_location"])
    })
}

values <- value(futures)
