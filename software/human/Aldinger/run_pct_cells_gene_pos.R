# ==============================================================================
# trying to answer 2 questions:
# 1. In Vladoiu et al, what % of cells are RL/UBC at each timepoint
# 2. In Aldinger et al, what % of cells in RL lineage are expressing a gene
# ==============================================================================

# set working directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
data.dir <- "/data/scRNAseq/"
# output.dir <- 

library(AnnotationHub)
library(Seurat)
library(tidyverse)

# ==============================================================================
# analyzing the Vladoiu data

# ------------------------------------------------------------------------------
# import the data


# ==============================================================================
# analyzing the Aldinger data

# ------------------------------------------------------------------------------
# import the data

expr.path <- paste0(data.dir, "CBL-dev/exprMatrix.tsv.gz")
meta.path <- paste0(data.dir, "CBL-dev/meta.tsv")

# import the expression matrix and the meta data
exprMatrix <- read.csv(expr.path, sep = "\t", row.names = 1)
metadata <- read.csv(meta.path, sep = "\t", row.names = 1)

# create Seurat object
human.srat <- CreateSeuratObject(exprMatrix, project = "CBL-dev", meta.data = metadata)

# convert gene expression to Z-score (scale and centre the data)
human.srat <- ScaleData(human.srat, features = rownames(human.srat), vars.to.regress = "CC.Difference")

expr <- human.srat[["RNA"]]@scale.data["EOMES", ]
timepoints <- human.srat@meta.data[, c("age", "Cluster")] %>% 
  mutate_if(is.character, as.factor)
eomes.expr <- cbind(timepoints, expr)

eomes.pct <- eomes.expr %>% 
  group_by(Cluster) %>% 
  group_map(~ PercentAbove(.x$expr, 0)) %>% 
  setNames(., levels(eomes.expr$Cluster)) %>% 
  unlist(.)
