library(tidyverse)
library(Seurat)
library(qs)

luo_srat <- qread("/.mounts/labs/pailab/private/llau/data/Luo_2022/Luo_RL_human_20240424_stand.qs")
Idents(luo_srat) <- "orig.ident"

luo_srat_subset <- subset(x = luo_srat, idents = c("FB13", "FB15"), invert = TRUE)


luo_srat_subset <- SCTransform(luo_srat_subset, vars.to.regress = "CC.Difference", variable.features.n = 5000, return.only.var.genes = FALSE) %>%
    RunPCA()


qsave(luo_srat_subset, "/.mounts/labs/pailab/private/llau/data/Luo_2022/Luo_RL_human_20240503_stand.qs")