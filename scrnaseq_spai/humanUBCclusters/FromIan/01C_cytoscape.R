# ==============================================================================
# Make EnrichmentMap in Cytoscape using RCy3.
#
# Note that this script needs to be run locally on a machine that has Cytoscape
# downloaded, NOT on the lab server or the HPC.
#
# SP: The current version is set to run on SP's laptop.
# ==============================================================================

library(argparse)
library(tidyverse)
library(RCy3)

rootDir <- "/Users/spai/Library/CloudStorage/GoogleDrive-spai@oicr.on.ca/Shared drives/Pai Lab/Projects/HumanMouseUBC"

# parse command line arguments
parser <- ArgumentParser()

parser$add_argument(
  # folder to g:Profiler results (MUST BE RELATIVE PATH)
  "--gost_dir",
  default = NULL,
  required = TRUE
)
parser$add_argument(
  # output directory
  "--out_dir",
  default = NULL,
  required = TRUE
)

dt <- format(Sys.Date(),"%y%m%d")

if (interactive()) {
  # for testing and troubleshooting
  arg_list <- parser$parse_args(c(
    "--gost_dir", sprintf("%s/HumanUBCclusters/PathwayEnrichment/gProfilerResults/20250401",rootDir),
    "--out_dir", sprintf("%s/HumanUBCclusters/PathwayEnrichment/EnrichmentMaps",rootDir)
  ))
} else {
  arg_list <- parser$parse_args()
}

arg_list$out_dir <- sprintf("%s/%s", arg_list$out_dir, dt)

message(sprintf("Saving files to %s", arg_list$out_dir))
if (!dir.exists(arg_list$out_dir)) {
  dir.create(arg_list$out_dir, recursive = FALSE)
}

# set full path to directories (necessary for Cytoscape)
for (d in c("gost_dir", "out_dir")) {
  arg_list[[d]] <- tools::file_path_as_absolute(arg_list[[d]])
}

# ------------------------------------------------------------------------------
# set parameters for EnrichmentMap (see
# https://enrichmentmap.readthedocs.io/en/latest/Parameters.html for more
# information)

# p-value threshold to filter genesets
pvalue_gprofiler_threshold <- 1.0
# q-value threshold to filter genesets
qvalue_gprofiler_threshold <- 0.05

# similarity threshold to filter genesets connections/edges
similarity_threshold <- "0.4"
# similarity metric to filter geneset connections/edges
similarity_metric <- "OVERLAP"

# GMT file
#gmt_file <- tools::file_path_as_absolute("./#Human_GOBP_AllPathways_no_GO_iea_August_08_2023_symbol_min10_max250.gmt")
gmt_file <- "/Users/spai/Documents/Pai_Lab/Projects/HumanMouseUBC/anno/Human_GOBP_AllPathways_no_GO_iea_August_08_2023_symbol_min10_max250.gmt"


# ------------------------------------------------------------------------------
# check Cytoscape connection

message("***Checking Cytoscape connection***")
cytoscapePing()
cytoscapeVersionInfo()

browser()

# ------------------------------------------------------------------------------
# create EnrichmentMap for each UBC cluster

ubc_subdir <- paste0("UBC_", c(1:5, 6, 7, 8))

walk(
  .x = ubc_subdir,
  .f = \(clust) {
    message(sprintf("***Generating EnrichmentMap for %s***", clust))

    em_command <- paste(
      'enrichmentmap build analysisType="generic"',
      'gmtFile=', gmt_file,
      'pvalue=', pvalue_gprofiler_threshold,
      'qvalue=', qvalue_gprofiler_threshold,
      'similaritycutoff=', similarity_threshold,
      'coefficients=', similarity_metric,
      'enrichmentsDataset1=', file.path(arg_list$gost_dir, clust, "enriched_pathways.gem.txt"),
      'filterByExpressions=false',
      sep = " "
    )

    # run em_command
    response <- commandsGET(em_command)

    # get suid of network
    if (grepl(pattern = "Failed", x = response)) {
      message(sprintf("***Failed to create network for %s", clust))
    } else {
      current_network_suid <- response
    }

    # rename network
    . <- renameNetwork(
      title = clust,
      network = as.numeric(current_network_suid)
    )

    # save screenshot
    fitContent()
    network_file <- file.path(arg_list$out_dir, clust)
    walk(
      .x = c("SVG", "PDF"),
      .f = \(x) {
        if(file.exists(paste0(network_file, x))){
          # cytoscape hangs waiting for user response if file already exists; remove it
          # first
          . <- file.remove(network_file)
        }
        . <- exportImage(filename = network_file, type = x)
      }
    )

    # save and close session
    . <- saveSession(filename = file.path(arg_list$out_dir, paste0(clust, ".cys")))
    closeSession(save.before.closing = FALSE)
  }
)

# ------------------------------------------------------------------------------

message("\n***SESSION INFO***\n")
print(sessionInfo())