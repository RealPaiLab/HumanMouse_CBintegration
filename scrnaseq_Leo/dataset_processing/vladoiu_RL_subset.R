library(Seurat)
library(tidyverse)
library(qs)

setwd("/u/llau/software/mb_scrnaseq/MB_scRNAseq")
source("scrnaseq_Leo/utilities/isolate_RL_lineage.R")

# Luo human FC
isolate_RL_lineage("/.mounts/labs/pailab/private/llau/data/Vladoiu_2019/Vladoiu_full_cerebellum_mouse_20240403_stand.qs", 
                    "/.mounts/labs/pailab/private/llau/data/Vladoiu_2019", 
                    "Vladoiu", 
                    "mouse")