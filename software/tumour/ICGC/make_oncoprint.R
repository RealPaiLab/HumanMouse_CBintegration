# ==============================================================================
# Make the actual Oncoprint from matrix.
# ==============================================================================

library(tidyverse)
library(ComplexHeatmap)

mut_file <- "/CBL_scRNAseq/results/tumour/ICGC/20231102/simple_somatic_mutation.PEME-CA.csv"

# ------------------------------------------------------------------------------

mut_mat <- read.csv(mut_file, row.names = 1) %>% as.matrix()

# colours
col <- c("snv" = "red", "indel <=200" = "blue")
alter_fun <- list(
  snv = alter_graphic(
    graphic = "rect",
    width = 0.9,
    height = 0.8,
    fill = col["snv"]
  ),
  `indel <=200` = alter_graphic(
    graphic = "rect",
    width = 0.9,
    height = 0.4,
    fill = col["indel <=200"]
  )
)

oncoPrint(
  mut_mat[1:100, ],
  alter_fun = alter_fun,
  col = col,
  column_title = "patient",
  row_title = "gene",
  show_row_names = FALSE,
  remove_empty_columns = TRUE,
  remove_empty_rows = TRUE,
  heatmap_legend_param = list(title = "alterations", at = c("snv", "indel <=200"))
)
