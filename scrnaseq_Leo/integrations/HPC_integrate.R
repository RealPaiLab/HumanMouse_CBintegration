# ==============================================================================
# Integrate datasets
# ==============================================================================

# ==============================================================================
# Create flags for arguments
# ==============================================================================

library(argparse)

# Parse command line arguments
parser <- ArgumentParser()

parser$add_argument(
  # Input dataset locations
  "-datasets",
  action = "store",
  nargs = "+",
  required = TRUE,
  help = 'List of names of datasets to be integrated in format {Author}_{Cells}_{Species}'
)

parser$add_argument(
  # Input intermediate integration file
  "-intermediate",
  action = "store",
  default = NULL,
  help = 'Location of intermediate integration file'
)

parser$add_argument(
  # Use all integration methods
  "--all_int",
  action = "store_true",
  default = FALSE,
  help = 'Integrate using all methods'
)

parser$add_argument(
  # Use CCA integration method
  "--CCA_int",
  action = "store_true",
  default = FALSE,
  help = 'Integrate using CCA method'
)

parser$add_argument(
  # Use Harmony integration method
  "--HARM_int",
  action = "store_true",
  default = FALSE,
  help = 'Integrate using Harmony method'
)

parser$add_argument(
  # Use RPCA integration method
  "--RPCA_int",
  action = "store_true",
  default = FALSE,
  help = 'Integrate using RPCA method'
)

parser$add_argument(
  # Do t-SNE plots
  "--tsne",
  action = "store_true",
  default = FALSE,
  help = 'Plot using t-SNE reduction'
)

parser$add_argument(
  # Run SCTransform before integration prep
  "--sctransform",
  action = "store_true",
  default = FALSE,
  help = 'Run SCTransform before integration prep'
)

parser$add_argument(
  # Using .rds files
  "--use_rds",
  action = "store_true",
  default = FALSE,
  help = 'Read/write .rds files for integration'
)

parser$add_argument(
  # Run pipeline without future()
  "--run_series",
  action = "store_true",
  default = FALSE,
  help = 'Run integrations in series'
)

arg_list <- parser$parse_args()

dataset_names <- arg_list$datasets
intermediate_file <- arg_list$intermediate
int_using_ALL <- arg_list$all_int
int_using_CCA <- arg_list$CCA_int
int_using_HARM <- arg_list$HARM_int
int_using_RPCA <- arg_list$RPCA_int
do_tsne <- arg_list$tsne
run_sctransform <- arg_list$sctransform
use_rds <- arg_list$use_rds
run_series <- arg_list$run_series


# ==============================================================================
# Loading in pipeline scripts and setting up directories
# ==============================================================================


# Setting working directory
setwd("/u/spai/software/MB_scRNAseq")

library(biomaRt)
library(tidyverse)
library(Seurat)
library(harmony)
library(RColorBrewer)
library(future)
library(qs)
#library(rliger)
library(SeuratWrappers)

# Importing integration pipeline
source("scrnaseq_Leo/pipeline_scripts/integration_prep.R")
source("scrnaseq_Leo/pipeline_scripts/cca_integration.R")
source("scrnaseq_Leo/pipeline_scripts/harmony_integration.R")
source("scrnaseq_Leo/pipeline_scripts/rpca_integration.R")


# Set data and output directories
integ_out_dir <- file.path("/.mounts/labs/pailab/private/llau", "results/integrated")
integ_date_dir <- file.path(integ_out_dir, format(Sys.Date(), "%Y%m%d"))
integ_CCA_dir <- file.path(integ_date_dir, "cca")
integ_HARM_dir <- file.path(integ_date_dir, "harmony")
integ_RPCA_dir <- file.path(integ_date_dir, "rpca")

if (!dir.exists(integ_date_dir) ) {
  dir.create(integ_date_dir)
}
if (!dir.exists(integ_CCA_dir) ) {
  dir.create(integ_CCA_dir)
}
if (!dir.exists(integ_HARM_dir) ) {
  dir.create(integ_HARM_dir)
}
if (!dir.exists(integ_RPCA_dir) ) {
  dir.create(integ_RPCA_dir)
}

# Resolutions to run FindClusters()
cluster_resolutions =  seq(0.2, 0.8, by = 0.2)

if (is.null(intermediate_file)) {
# ==============================================================================
# Import data and fill in parameters
# ==============================================================================

  datasets <- read.csv("scrnaseq_Leo/integrations/dataset_list.csv")
  # Read each RDS file and store its content in the list
  all_datasets <- list()
  for(i in seq_along(dataset_names)){
    file <- subset(datasets, dataset_name == dataset_names[i])
    # Checking to see if user input is correct
    message("Loading in: ", dataset_names[i])
    if(nrow(file) == 0){
      stop("Dataset does not exist. Please check integrations/dataset_list.csv to see if specified dataset is on the list.")
    }
    if(use_rds){
      message("Loading in: ", file$dataset_location_rds)
      all_datasets[[i]] <- readRDS(file$dataset_location_rds)
    } else {
      message("Loading in: ", file$dataset_location_qs)
      all_datasets[[i]] <- qread(file$dataset_location_qs)
    }
  }

  # ==============================================================================
  # Prepping to integrate the datasets
  # ==============================================================================

  # A list containing the integrated datasets and integration features
  integration_vars <- integration_prep(all_datasets, run_sctransform)

  # Saving intermediate integration file
  if(use_rds){
    print("Saving intermediate integration file as .rds file")
    file_name = paste(format(Sys.Date(), "%Y%m%d"), "intermediate_integ.rds", sep = "_")
    saveRDS(object = integration_vars, file = file.path(integ_date_dir, file_name))
  } else {
    print("Saving intermediate integration file as .qs file")
    file_name = paste(format(Sys.Date(), "%Y%m%d"), "intermediate_integ.qs", sep = "_")
    qsave(integration_vars, file.path(integ_date_dir, file_name))
  }

  # Memory management
  rm(all_datasets)
  gc()
} else {
  # Using premade intermediate integration file
  print("Reading in intermediate integration file")
  if(use_rds){
    integration_vars <- readRDS(file.path(intermediate_file))
  } else {
    integration_vars <- qread(file.path(intermediate_file))
  }
}

options(future.globals.maxSize = 200 * 1024^3)
message("Available cores: ", as.integer(availableCores()))
#options(future.globals.onReference = "error")

if(run_series){
  # ==============================================================================
  # CCA Integration
  # ==============================================================================

  if(int_using_ALL | int_using_CCA){
    print("Beginning CCA integration")
    cca_integration( integ_list = integration_vars, resolutions = cluster_resolutions, 
              do_tsne = do_tsne, out_directory = integ_CCA_dir, use_rds = use_rds)
    gc()
  }

  # ==============================================================================
  # Harmony Integration
  # ==============================================================================

  if(int_using_ALL | int_using_HARM){
    print("Beginning Harmony integration")
    harmony_integration( integ_list = integration_vars, resolutions = cluster_resolutions, 
              do_tsne = do_tsne, out_directory = integ_HARM_dir, use_rds = use_rds)
    gc()
  }
 
  # ==============================================================================
  # RPCA Integration
  # ==============================================================================

  if(int_using_ALL | int_using_RPCA){
    print("Beginning RPCA integration")
    rpca_integration( integ_list = integration_vars, resolutions = cluster_resolutions, 
            do_tsne = do_tsne, out_directory = integ_RPCA_dir, use_rds = use_rds)
    gc()
  }

  message("Integration is complete")
} else {
  # Enable parallelization
  plan("multicore", workers = 3)

  cca %<-% {
      # ==============================================================================
      # CCA Integration
      # ==============================================================================

      if(int_using_ALL | int_using_CCA){
        # Creating log file
        log_file <- file.path(integ_CCA_dir, "CCA_log_file.log")
        file_conn <- file(log_file, open = "a")
        sink(file_conn, type = "output")
        sink(file_conn, type = "message", append = TRUE)

        print("Beginning CCA integration")
        cca_integration( integ_list = integration_vars, resolutions = cluster_resolutions, 
                  do_tsne = do_tsne, out_directory = integ_CCA_dir, log_file = file_conn, use_rds = use_rds)

        sink(type = "output")
        sink(type = "message")
      }
      gc()
  }

  harm %<-% {
      # ==============================================================================
      # Harmony Integration
      # ==============================================================================

      if(int_using_ALL | int_using_HARM){
        # Creating log file
        log_file <- file.path(integ_HARM_dir, "Harmony_log_file.log")
        file_conn <- file(log_file, open = "a")
        sink(file_conn, type = "output")
        sink(file_conn, type = "message", append = TRUE)

        print("Beginning Harmony integration")
        harmony_integration( integ_list = integration_vars, resolutions = cluster_resolutions, 
                  do_tsne = do_tsne, out_directory = integ_HARM_dir, log_file = file_conn, use_rds = use_rds)

        sink(type = "output")
        sink(type = "message")
      }
      gc()
  }
      
  rpca %<-% {
      # ==============================================================================
      # RPCA Integration
      # ==============================================================================

      if(int_using_ALL | int_using_RPCA){
        # Creating log file
        log_file <- file.path(integ_RPCA_dir, "RPCA_log_file.log")
        file_conn <- file(log_file, open = "a")
        sink(file_conn, type = "output")
        sink(file_conn, type = "message", append = TRUE)

        print("Beginning RPCA integration")
        rpca_integration( integ_list = integration_vars, resolutions = cluster_resolutions, 
                do_tsne = do_tsne, out_directory = integ_RPCA_dir, log_file = file_conn, use_rds = use_rds)

        sink(type = "output")
        sink(type = "message")
      }
      gc()
  }

  cca_resolve <- futureOf(cca)
  harm_resolve <- futureOf(harm)
  rpca_resolve <- futureOf(rpca)

  time <- 0
  while (!resolved(cca_resolve) | !resolved(harm_resolve) | !resolved(rpca_resolve)) {
    message("Integration under way. Currently at " , (time/60), " minutes")
    Sys.sleep(60)
    time <- time + 60
  }

  message("Integration is complete")
}

print(sessionInfo())


