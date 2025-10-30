# check expression of PRDM6, GFI1, MAX, E2F6, TFDP1, FOXP2 in Aldinger data

setwd(dirname(rstudioapi::getSourceEditorContext()$path))

library(ggplot2)
library(Seurat)

genes <- c("PRDM6", "GFI1", "MAX", "E2F6", "TFDP1", "FOXP2")

srat <- readRDS("/isilon/CBL_scRNAseq-archived/data/human/Aldinger/glutamatergic_dev_Liam.RDS")

# normalize RNA assay (results saved to `srat[["RNA"]]@data` slot)
srat <- NormalizeData(srat, assay = "RNA")

.plt <- FeaturePlot(
  srat,
  features = genes,
  split.by = "new_cell_type",
  keep.scale = "all",
  # min.cutoff = "q10",
  # max.cutoff = "q90",
  order = TRUE
)
ggsave(
  "gene_expression.png",
  plot = .plt,
  width = 20,
  height = 20,
  units = "in"
)

# normalize counts then 

avg_exp <- AverageExpression(
  srat,
  features = genes,
  group.by = "new_cell_type"
)
write.csv(
  avg_exp$RNA,
  file = "average_expression.csv",
)
