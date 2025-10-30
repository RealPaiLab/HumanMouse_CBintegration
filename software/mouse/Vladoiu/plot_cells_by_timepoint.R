# ==============================================================================
# Take the Vladoiu data and plot the number of cells from each timepoint
# ==============================================================================

# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
# setwd("/CBL_scRNAseq/software/")

library(tidyverse)
library(magrittr)
library(Seurat)

# set data and output directories
root_dir <- "CBL_scRNAseq"
data_dir <- file.path("/isilon", root_dir, "data/mouse/Vladoiu")
out_dir <- file.path("", root_dir, "results/mouse/Vladoiu")
date_dir <- file.path(out_dir, format(Sys.Date(), "%Y%m%d"))

# create output directory if it doesn't exist
if (!dir.exists(date_dir)) {
  dir.create(date_dir)
}

# load Seurat object from RDS file
srat <- readRDS(file.path(out_dir, "merged_seurat.rds"))

# get number of cells in each timepoint and convert into a data frame
timepoints <- table(srat$orig.ident) %>% 
  data.frame(.) %>% 
  set_colnames(c("timepoint", "num_cells")) %>% 
  mutate(timepoint = str_remove(timepoint, "Vladoiu-"), 
         timepoint = fct_relevel(timepoint, "P14", after = Inf)) %>% 
  arrange(timepoint)

# make and save plot
plt <- ggplot(timepoints, aes(x = timepoint, y = num_cells, fill = timepoint)) + 
  geom_col(width = 0.6) + 
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) + 
  theme_classic() +
  theme(axis.text = element_text(colour = "black"))
ggsave(filename = "cells_by_timepoint.pdf", plot = plt, path = date_dir, 
       width = 6, height = 4, units = "in")
