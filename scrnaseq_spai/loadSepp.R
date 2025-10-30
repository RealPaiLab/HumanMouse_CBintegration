
library(Seurat)

inFile <- "/home/rstudio/isilon/private/llau/data/Sepp_2024/Sepp_full_cerebellum_human_20240509_stand.rds"

srat <- readRDS(inFile)
mdata <- srat@meta.data

outFile <- "/home/rstudio/isilon/private/projects/HumanMouseUBC/Sepp2024_metadata.txt"
write.table(mdata,file=outFile,sep="\t",col=TRUE,row=TRUE,quote=FALSE)
