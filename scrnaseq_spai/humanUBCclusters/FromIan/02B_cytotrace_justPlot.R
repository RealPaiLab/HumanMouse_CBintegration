rm(list=ls())
library(tidyverse)
library(Seurat)

# set python for CytoTRACE
#reticulate::use_condaenv("monocle3")
#library(CytoTRACE)

 srat_file <- "/home/rstudio/isilon/private/projects/HumanMouseUBC/UBCclusters/fromQuang/UBC.Harmony.RDS"

inDir <- "/home/rstudio/isilon/private/projects/HumanMouseUBC/integrated_human_ubc/CytoTRACE/250424/with_dataset_correction"

inFile <- sprintf("%s/cytotrace_output.rds", inDir)

colourFile <- "../UBCcolours.txt"
clrs <- read.delim(colourFile,header=TRUE,sep="\t")
   

cat("Reading seurat file\n")
t0 <- Sys.time()
srat <- readRDS(srat_file)
print(Sys.time() - t0)

cat("Reading CytoTRACE output\n")
t0 <- Sys.time()
ct_res <- readRDS(inFile)
print(Sys.time() - t0)

dt <- format(Sys.Date(),"%y%m%d")
outDir <- sprintf("%s/%s",inDir,dt)
if (!dir.exists(outDir)) {
  dir.create(outDir, recursive = FALSE)
}

srat$cytotrace <- ct_res$CytoTRACE[rownames(srat[[]])]    

md <- srat[[]]
md$SCT_snn_res.0.5 <- factor(md$SCT_snn_res.0.5,levels=clrs$Cluster)

md<- md[-which(md$SCT_snn_res.0.5 %in% 6),]
clrs <- clrs[-which(clrs$Cluster %in% 6),]
# filename for plot
filename <- sprintf("vln_%s.pdf", dt)
p <- ggplot(md, aes(x = SCT_snn_res.0.5, y = cytotrace)) + 
  geom_violin(aes(fill = SCT_snn_res.0.5)) + 
  geom_jitter(width = 0.25, size = 0.1) + 
  scale_fill_manual(values = clrs$Colour) +
  labs(x = "Clusters", y = "CytoTRACE score\n(higher = more stem-like)") + 
  theme_classic(base_size=21) + 
  theme(
    axis.text = element_text(colour = "black"),
    axis.ticks = element_line(colour = "black"),
    legend.position = "none"
  )
# save plot
ggsave(p, file=sprintf("%s/%s",outDir,filename),
  width = 8,
  height = 4,
  units = "in"
)
