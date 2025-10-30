library(tidyverse)
library(Seurat)
library(qs)

aldinger_srat <- qread("/.mounts/labs/pailab/private/llau/data/Aldinger_2021/Aldinger_RL_human_20240415_stand.qs")

Idents(aldinger_srat) <- "age"
aldinger_srat_subset <- subset(x = aldinger_srat, idents = c("9 PCW", "10 PCW"), invert = TRUE)

aldinger_srat_subset <- SCTransform(aldinger_srat_subset, vars.to.regress = "CC.Difference", variable.features.n = 5000, return.only.var.genes = FALSE) %>%
    RunPCA()

qsave(aldinger_srat_subset, "/.mounts/labs/pailab/private/llau/data/Aldinger_2021/Aldinger_RL_human_20240514_stand.qs")