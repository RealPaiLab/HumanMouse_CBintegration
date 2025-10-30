# create ComplexHeatmap of Erickson clusters using top UBC cluster genes

library(ComplexHeatmap)
library(ggplot2)
library(reshape2)
library(ggpubr)

mb <- "/home/rstudio/isilon/private/projects/MB_multiome/output/QC/20250206/fastmnn/mb_fastmnn.qs"
outRoot <- "/home/rstudio/isilon/private/projects/HumanMouseUBC/EricksonHeatmap"

outDir <- sprintf("%s/%s", outRoot, format(Sys.Date(), "%y%m%d"))
if (!dir.exists(outDir)) {
    dir.create(outDir, recursive = FALSE)
}

clusterFile <- "/home/rstudio/isilon/private/projects/MB_multiome/output/clustering/250502/clustering_assignments.txt"

cl <- read.table(clusterFile, header=TRUE,sep="\t")
 if (all.equal(rownames(tumour_srat[[]]),rownames(cl))!=TRUE){
     cat("Row names of metadata and clustering do not match. Fix this. \n")
     browser()
 }
 tumour_srat[["snn_res.1"]] <- factor(cl[,"snn_res.1"])
 tumour_srat$seurat_clusters <- factor(cl[,"snn_res.1"])


cat("Loading Erickson MB data - takes 13 min on Pai lab server\n")
t0 <- Sys.time()
srat <- qs::qread(mb)
print(Sys.time() - t0)

# Compute AverageExpression of all genes for the Erickson dataset.
cat("Computing average expression of all genes\n")

x <- c("MKI67","TOP2A","E2F7","HMGA2","RARB","E2F3",
       "PBX3","PAX2","SLC6A1","ITPR1","SORCS3",
       "SMAD9","PLAGL1","GREM2","LINC00486","CTNT2",
       "EYS","ICQJ-SCHIP1",
       "PEX5L",
       "RFX3","TMEM163", "MGAT4C")

# take intersection of UBC genes and top genes from the Erickson dataset
y <- intersect(x,rownames(srat))

xpr <- AggregateExpression(
  tumour_srat,
  assays = "SCT",
  features = y,
  group.by = "snn_res.1",
  slot = "data"
)[[1]]

xpr <- t(scale(t(xpr)))

# heatmap colours
col_fun <- circlize::colorRamp2(
  breaks = c(-2,0,2),
  colors = c("#FF00FF","#000000", "#FFFF00")
)

# plot heatmap
cat("Plotting heatmap\n")
hm <- Heatmap(
  xpr,
  name="scaled\nexpression",
  col = col_fun,
)

pdf(
  file = file.path(outDir, "heatmap.pdf"),
  width = 14,
  height = 20
)
draw(hm, annotation_legend_list = lgd_list, merge_legend = TRUE)
dev.off()