rm(list=ls())

library(Seurat)

inDir <- "/home/rstudio/isilon/src/neurodev-genomics/SpatialRNA/GSE165657/GSM5858166_10xspatial_GW17_05"
outDir <- "/home/rstudio/isilon/private/projects/HumanMouseUBC/ZhongSpatial"

dt <- format(Sys.Date(),"%y%m%d")
outDir <- sprintf("%s/%s", outDir, dt)
if (!dir.exists(outDir)) {
  dir.create(outDir, recursive = T)
}

t0 <- Sys.time()
obj <- Read10X(data.dir = inDir)
print(Sys.time()-t0)

# create seurat object from 10X data
srat <- CreateSeuratObject(
    counts = obj,
    project = "ZhongSpatial",
    min.cells = 3,
    min.features = 200
)

srat[["percent.mt"]] <- PercentageFeatureSet(srat, pattern = "^MT-")
