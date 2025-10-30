library(tidyverse)
library(Seurat)
library(qs)

sepp_srat <- qread("/.mounts/labs/pailab/private/llau/data/Sepp_2024/Sepp_RL_human_20240415_stand.qs")

Idents(sepp_srat) <- "Stage"
sepp_srat_subset <- subset(x = sepp_srat, idents = c("7 wpc", "8 wpc", "9 wpc"), invert = TRUE)

sepp_srat_subset <- SCTransform(sepp_srat_subset, vars.to.regress = "CC.Difference", variable.features.n = 5000, return.only.var.genes = FALSE) %>%
    RunPCA()

qsave(sepp_srat_subset, "/.mounts/labs/pailab/private/llau/data/Sepp_2024/Sepp_RL_human_20240514_stand.qs")