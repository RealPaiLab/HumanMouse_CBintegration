library(Seurat)
library(ggplot2)
library(pheatmap)

# devtools::install_github("dmcable/spacexr", build_vignettes = FALSE)
library(spacexr)


outDir <- "/.mounts/labs/pailab/private/xsun/output/hb_visium_deconv"

## data prep ##
visium <- readRDS("/.mounts/labs/pailab/private/projects/FetalHindbrain/Aldinger_Visium_2021Data/obj.list.rds")

# check slices
#SpatialFeaturePlot(visium[[1]], features = "EOMES", alpha = 0.5)
#SpatialFeaturePlot(visium[[2]], features = "EOMES", alpha = 0.5)
#SpatialFeaturePlot(visium[[3]], features = "EOMES", alpha = 0.5)
#SpatialFeaturePlot(visium[[4]], features = "EOMES", alpha = 0.5)

visium[[1]]@meta.data$orig.ident <- "sample1_slice1"
visium[[2]]@meta.data$orig.ident <- "sample1_slice2"
visium[[3]]@meta.data$orig.ident <- "sample2_slice1"
visium[[4]]@meta.data$orig.ident <- "sample2_slice2"


# Leo & Ian's UBC annotation
ubc_annotation_file <- "/.mounts/labs/pailab/private/xsun/tmp/ubc_subset_metadata.csv"
ubc_annotation <- read.csv(ubc_annotation_file, stringsAsFactors = F, row.names = 1)
ubc_annotation$ref_cell_type <- paste("UBC", ubc_annotation$snn_res.0.3, sep = "_")

# Aldinger full cerebellum #
# use @ to get/set meta because the object doesn't fit Seurat v5
aldinger <- readRDS("/.mounts/labs/pailab/src/neurodev-genomics/scRNAseq/Aldinger_2021/seurat.rds")

# transfer ubc cluster info
aldinger@meta.data$ref_cell_type <- gsub("H-", "", aldinger@meta.data$fig_cell_type)
tmp <- ubc_annotation$ref_cell_type
names(tmp) <- rownames(ubc_annotation)
tmp <- tmp[intersect(rownames(aldinger@meta.data), names(tmp))]
aldinger@meta.data[names(tmp),]$ref_cell_type <- tmp
print(table(aldinger@meta.data$ref_cell_type))

hm <- pheatmap(pmin(table(aldinger@meta.data$ref_cell_type, aldinger@meta.data$fig_cell_type), 100),
               scale = "none",
               fontsize = 14,
               fontsize_row = 14,
               fontsize_col = 14,
               color = colorRampPalette(c("white", "firebrick3"))(50),
               main = "Heatmap comparing labels (n capped at 100)", 
               cluster_rows = F, cluster_cols = F
)
ggsave(sprintf("%s/test_visium_4/aldingerUBC1_preCellType_vs_cellType_heatmap.png", outDir),
       hm, 
       dpi = 600, width = 20, height = 20
)

aldinger@meta.data$ref_cell_type[aldinger@meta.data$ref_cell_type == "eCN/UBC"] <- NA
aldinger@meta.data$ref_cell_type <- gsub("/", "_", aldinger@meta.data$ref_cell_type)

# Sepp full cerebellum #
sepp <- readRDS("/.mounts/labs/pailab/private/xsun/tmp/sepp2024_human_sct.rds")
sepp <- subset(sepp, cells = colnames(sepp)[!is.na(sepp$cell_type)])

# transfer UBC 0-5 cluster names to sepp
sepp$ref_cell_type <- sepp$cell_type
sepp <- AddMetaData(sepp, ubc_annotation[, c("snn_res.0.3", "ref_cell_type")])
print(table(sepp$ref_cell_type))

hm <- pheatmap(pmin(table(sepp$ref_cell_type, sepp$precisest_label), 100),
               scale = "none",
               fontsize = 14,
               fontsize_row = 14,
               fontsize_col = 14,
               color = colorRampPalette(c("white", "firebrick3"))(50),
               main = "Heatmap comparing labels (n capped at 100)", 
               cluster_rows = F, cluster_cols = F
)
ggsave(sprintf("%s/test_visium_4/seppUBC1_preciseCellType_vs_cellType_heatmap.png", outDir),
       hm, 
       dpi = 600, width = 20, height = 20
)

# remove leftover sepp_UBC 
sepp <- subset(sepp, subset = ref_cell_type != "UBC")
#sepp$ref_cell_type[sepp$ref_cell_type == "UBC"] <- NA
print(table(sepp$ref_cell_type))

# format cell type name as RCTD doesn't allow "/"
sepp$ref_cell_type <- gsub("/", "_", sepp$ref_cell_type)

# combine all UBC cluster into one
sepp$ref_cell_type_combineUBC <- sepp$ref_cell_type
sepp$ref_cell_type_combineUBC[startsWith(sepp$ref_cell_type_combineUBC, "UBC_")] <- "UBC"

# further separate GC to include GCP; 
sepp$ref_cell_type_precise <- sepp$ref_cell_type
sepp$ref_cell_type_precise[sepp$precisest_label == "GCP" & (! startsWith(
  sepp$ref_cell_type_precise, 
  "UBC_")
  )] <- "GCP"
print(table(sepp$ref_cell_type_precise))

# also test precise with gc/ubc cells removed
sepp$ref_cell_type_precise_noGCUBC <- sepp$ref_cell_type_precise
sepp$ref_cell_type_precise_noGCUBC[sepp$ref_cell_type_precise_noGCUBC == "GC_UBC"] <- NA

## testing RCTD reference settings ##

#' deconvolute hindbrain visium against scRNA ref using RCTD
#' @param visium (SeuratObject) visium spatial object
#' @param ref (SeuratObject) scRNA reference
#' @param ref_cell_type (character) Column to be used as reference cell type
#' @param max_cores (numeric) Max number of cores for RCTD
#' @return (SeuratObject) The original visium object with deconvoluted cell type weigths
deconv_hb_visium_RCTD <- function(visium, 
                                  ref, ref_cell_type, 
                                  max_cores = 8) {
  # extract information to pass to the RCTD Reference function
  counts <- ref@assays$RNA$counts
  cluster <- as.factor(ref@meta.data[[ref_cell_type]])
  names(cluster) <- rownames(ref@meta.data)
  nUMI <- ref@meta.data$nCount_RNA
  names(nUMI) <- rownames(ref@meta.data)
  reference <- Reference(counts, cluster, nUMI)
  
  # set up query with the RCTD function SpatialRNA
  counts <- visium[["Spatial"]]$counts
  coords <- GetTissueCoordinates(visium)
  colnames(coords) <- c("x", "y")
  coords[is.na(colnames(coords))] <- NULL
  query <- SpatialRNA(coords, counts, colSums(counts))
  
  # RCTD
  RCTD <- create.RCTD(query, reference, max_cores = max_cores)
  RCTD <- run.RCTD(RCTD, doublet_mode = "full")
  
  # assign cell type weights
  weights <- RCTD@results$weights
  norm_weights <- normalize_weights(weights)
  
  visium <- AddMetaData(visium, metadata = norm_weights)
  for (col in colnames(weights)) {
    visium@meta.data[[col]][is.na(visium@meta.data[[col]])] <- 0
  }
  
  return(visium)
}

# Test to find appropriate ref settings #
if (F) {
  # ref_cell_type #
  test_clusterUBC <- deconv_hb_visium_RCTD(
    visium = visium[[4]], 
    ref = sepp, ref_cell_type = "ref_cell_type", 
    max_cores = 10
  )
  
  clusterUBC_SpatialPlot <- SpatialPlot(test_clusterUBC, 
                                        features = unique(sepp$ref_cell_type), 
                                        alpha = 0.7, keep.scale = "all")
  
  ggsave(sprintf("%s/test_visium_4/clusterUBC_SpatialPlot.png", outDir), 
         clusterUBC_SpatialPlot, width = 40, height = 35, dpi = 600)
  
  
  # ref_cell_type_combineUBC #
  # first deconv with broader UBC 
  test_combinedUBC <- deconv_hb_visium_RCTD(
    visium = visium[[4]], 
    ref = sepp, ref_cell_type = "ref_cell_type_combineUBC", 
    max_cores = 10
  )
  
  
  combinedUBC_SpatialPlot <- SpatialPlot(test_combinedUBC, 
                                         features = unique(sepp$ref_cell_type_combineUBC), 
                                         alpha = 0.7, keep.scale = "all")
  
  ggsave(sprintf("%s/test_visium_4/combinedUBC_SpatialPlot.png", outDir), 
         combinedUBC_SpatialPlot, width = 40, height = 35, dpi = 600)
  
  
  # precise label sepp #
  test_clusterUBC_precise <- deconv_hb_visium_RCTD(
    visium = visium[[4]], 
    ref = sepp, ref_cell_type = "ref_cell_type_precise", 
    max_cores = 16
  )
  
  clusterUBC_precise_SpatialPlot <- SpatialPlot(test_clusterUBC_precise, 
                                                features = unique(sepp$ref_cell_type_precise), 
                                                alpha = 0.7, keep.scale = "all")
  
  ggsave(sprintf("%s/test_visium_4/clusterUBC_precise_SpatialPlot.png", outDir), 
         clusterUBC_precise_SpatialPlot, width = 40, height = 35, dpi = 600)
  
  # precise label sepp remove ubc/gcp #
  test_clusterUBC_preciseWithoutGCUBC <- deconv_hb_visium_RCTD(
    visium = visium[[4]], 
    ref = sepp, ref_cell_type = "ref_cell_type_precise_noGCUBC", 
    max_cores = 20
  )
  
  clusterUBC_precise_SpatialPlot <- SpatialPlot(
    test_clusterUBC_preciseWithoutGCUBC, 
    features = unique(na.omit(sepp$ref_cell_type_precise_noGCUBC)), 
    alpha = 0.7, keep.scale = "all")
  
  ggsave(sprintf("%s/test_visium_4/clusterUBC_preciseWithoutGCUBC_SpatialPlot.png", outDir), 
         clusterUBC_precise_SpatialPlot, width = 40, height = 35, dpi = 600)
  
  
  # aldinger with UBC clusters #
  test_aldinger <- deconv_hb_visium_RCTD(
    visium = visium[[4]],
    ref = aldinger, ref_cell_type = "ref_cell_type",
    max_cores = 16
  )
  
  aldinger_SpatialPlot <- SpatialPlot(
    test_aldinger, 
    features = unique(na.omit(aldinger@meta.data$ref_cell_type)), 
    alpha = 0.7, keep.scale = "all")
  
  ggsave(sprintf("%s/test_visium_4/aldinger_SpatialPlot.png", outDir), 
         aldinger_SpatialPlot, width = 40, height = 35, dpi = 600)
}


## Generate for all sections ##
# save ref file for future usage #
dt <- format(Sys.Date(),"%y%m%d")

saveRDS(aldinger, sprintf("%s/ref_rds/aldinger_%s.rds", outDir, dt))
saveRDS(sepp, sprintf("%s/ref_rds/sepp_%s.rds", outDir, dt))

# Aldinger ref #
for (i in 1:4) {
  message(sprintf("--- Deconvoluting visium %i against Aldinger ref", i))
  
  tmp <- deconv_hb_visium_RCTD(
    visium = visium[[i]],
    ref = aldinger, 
    ref_cell_type = "ref_cell_type",
    max_cores = 16
  )
  
  .plot <- SpatialPlot(
    tmp, 
    features = unique(na.omit(aldinger@meta.data$ref_cell_type)), 
    alpha = 0.7, keep.scale = "all")
  
  saveRDS(tmp, 
          sprintf("%s/results/aldinger_ref/visium_%i_against_aldinger_%s.rds", 
                  outDir, i, dt)
          )
  
  ggsave(
    sprintf(
      "%s/results/aldinger_ref/visium_%i_against_aldinger_SpatialPlot_%s.png", 
      outDir, i, dt), 
    .plot, 
    width = 40, height = 35, dpi = 600)
  
  gc()
}

# Sepp ref #
for (i in 1:4) {
  message(sprintf("--- Deconvoluting visium %i against Sepp ref", i))
  
  tmp <- deconv_hb_visium_RCTD(
    visium = visium[[i]],
    ref = sepp, 
    ref_cell_type = "ref_cell_type_precise",
    max_cores = 16
  )
  
  .plot <- SpatialPlot(
    tmp, 
    features = unique(sepp$ref_cell_type_precise), 
    alpha = 0.7, keep.scale = "all")
  
  saveRDS(tmp, 
          sprintf("%s/results/sepp_ref/visium_%i_against_sepp_%s.rds", 
                  outDir, i, dt)
  )
  
  ggsave(
    sprintf(
      "%s/results/sepp_ref/visium_%i_against_sepp_SpatialPlot_%s.png", 
      outDir, i, dt), 
    .plot, 
    width = 40, height = 35, dpi = 600)
  
  gc()
}



