# ==============================================================================
# function to load Vladoiu count matrices
# ==============================================================================

# load Vladoiu count matrix as Seurat object
load_mtx <- function(
  timepoint, 
  data_dir, 
  mtx = "matrix.mtx.gz", 
  cells = "barcodes.tsv.gz", 
  features = "features.tsv.gz"
) {
  # read in the count matrix
  count_files <- file.path(data_dir, timepoint, c(mtx, cells, features))
  message(sprintf("* Reading in the count matrix for %s", timepoint))
  count_data <- ReadMtx(mtx = count_files[1],
                        cells = count_files[2],
                        features = count_files[3],
                        strip.suffix = TRUE)
  
  # create Seurat object
  message(sprintf("* Creating Seurat object for %s", timepoint))
  srat <- CreateSeuratObject(counts = count_data,
                             project = paste0("Vladoiu-", timepoint),
                             min.cells = 3,
                             min.features = 200)
   
  # rename cell barcodes to include timepoint
  new_names <- paste0(timepoint, "_", colnames(srat))
  srat <- RenameCells(srat, new.names = new_names)
  
  # copy barcodes to metadata
  srat@meta.data$barcode <- rownames(srat@meta.data)
  
  return(srat)
}

# ==============================================================================
# functions to filter cells
# ==============================================================================

# filter cells
filter_cells <- function (
  srat, 
  deviation = sd, 
  deviation_cutoff = 5
) {
  # calculate percentage of mitochondrial reads
  srat[["percent_mt"]] <- PercentageFeatureSet(srat, pattern = "^mt-")
  
  # calculate the deviation cutoff for total cell count (nCount_RNA) and percent_mt
  count_cut <- deviation(srat[["nCount_RNA"]][, 1]) * deviation_cutoff
  mt_cut <- deviation(srat[["percent_mt"]][, 1]) * deviation_cutoff
  
  srat <- subset(srat, 
                 subset = nCount_RNA < count_cut & percent_mt < mt_cut)
  
  return(srat)
}
