rm(list=ls())

library(tidyverse)
library(patchwork)
library(Seurat)
library(pals)
library(readxl)
source("FromIan/utils.R")

srat_qs <- get_srat_paths()
ubc_srat <- load_srat(srat_qs["ubc"]) %>%
  pluck("ubc")
##message("done")
my_pals <- get_custom_pals()


humanUBC <- subset(ubc_srat, subset = species == "human")
mdata <- humanUBC@meta.data
mdata$cell_id <- rownames(mdata)

out_dir <- "/home/rstudio/isilon/private/projects/HumanMouseUBC/UBC1_marker_analysis"
write.table(mdata,file=sprintf("%s/humanUBC_CellID.txt",out_dir),
    sep="\t",col=T,row=F,quote=F)