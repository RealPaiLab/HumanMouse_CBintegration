# ==============================================================================
# Get common TFs across cell types.
# ==============================================================================

# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(argparse)
library(tidyverse)
library(ggVennDiagram)
library(patchwork)

# parse command line arguments
parser <- ArgumentParser()

parser$add_argument(
  # path to CSV of RSS scores
  "--rss_csv",
  default = NULL,
  required = TRUE
)
parser$add_argument(
  # number of regulons to compare
  "--num_regulons",
  default = NULL,
  required = TRUE,
  type = "integer"
)
parser$add_argument(
  # output directory
  "--out_dir",
  default = NULL,
  required = TRUE
)

if (interactive()) {
  # for testing and troubleshooting
  args <- parser$parse_args(c(
    "--rss_csv", "/CBL_scRNAseq/results/pyscenic/20230320/aldinger_RL.rss.csv",
    "--num_regulons", 8,
    "--out_dir", "/CBL_scRNAseq/results/pyscenic/20230413/integ_cell_types/"
  ))
} else {
  args <- parser$parse_args()
}

message(sprintf("Saving files to %s", args$out_dir))
if (!dir.exists(args$out_dir)) {
  dir.create(args$out_dir, recursive = TRUE)
}

# ------------------------------------------------------------------------------
# load top pySCENIC regulons using the regulon specificity scores

# import data
rss <- read.csv(args$rss_csv, row.names = 1)

cell_types <- rownames(rss)

# top regulons for each cluster/cell type
top_regulons <- lapply(
  X = as.list(cell_types),
  FUN = function(x, rss, num_regulons) {
    # get the RSS scores
    scores <- data.matrix(rss)[x, ]
    
    # sort scores then get the TF associated with the score
    top_regs <- scores %>% 
      sort(decreasing = TRUE) %>% 
      names() %>% 
      head(n = num_regulons)
    
    return(top_regs)
  },
  rss = rss, # RSS scores
  num_regulons = args$num_regulons # top n regulons
)

names(top_regulons) <- cell_types

# ------------------------------------------------------------------------------
# Venn diagram for common TFs in cerebellar development and mutated genes in MB.

cell_type_pairs <- t(combn(cell_types, m = 2))

venn_plts <- pmap(
  .l = list(cell_type_pairs[,1], cell_type_pairs[,2]),
  .f = function(x, y, top_regulons) {
    # subset top TFs for this specific cell type
    tfs <- top_regulons[c(x, y)]
    
    # combine TFs and mutated genes into list and prep data for plotting
    venn_data <- list(tfs[[1]], tfs[[2]]) %>% 
      `names<-`(c(x, y)) %>% 
      Venn() %>% 
      process_data()
    
    # get intersection of TFs and mutated genes
    common_genes <- dplyr::intersect(tfs[[1]], tfs[[2]])
    
    # make the plot
    .plt <- ggplot() + 
      
      # the set areas
      geom_sf(aes(fill = count), data = venn_region(venn_data)) + 
      
      # the set edges
      geom_sf(color = "black", size = 2, data = venn_setedge(venn_data), show.legend = FALSE) + 
      
      # name of the sets
      geom_sf_text(
        aes(label = name),
        data = venn_setlabel(venn_data),
        size = 6,
        fontface = "bold"
      ) +
      
      # counts
      geom_sf_text(aes(label = count), data = venn_region(venn_data), size = 5) + 
      
      # label intersecting genes
      geom_segment(
        aes(x = 500, y = 400, xend = 500, yend = 200),
        color = "red",
        linewidth = 1
      ) + 
      geom_text(
        aes(x = 500, y = 180, label = paste(common_genes, collapse = "\n")),
        data = venn_region(venn_data),
        size = 2,
        vjust = 1
      ) + 
      
      # additional customizations
      scale_x_continuous(expand = expansion(mult = 0.5)) + 
      scale_y_continuous(expand = expansion(mult = c(0.5, 0.1))) + 
      scale_fill_distiller(palette = "Blues", direction = 1, guide = "none") + 
      theme_void()
    
    return(.plt)
  },
  top_regulons = top_regulons
)

# plt_dim <- length(cell_types) - 1
# plt_coords <- expand_grid(1:6, 2:7)
# areas <- area(t = plt_coords[,1], l = plt_coords)
# for (row in 1:nrow(plt_coords)) {
#   # paste(plt_coords[[row, 1]], plt_coords[[row, 2]]) %>% print()
#   a <- area(t = plt_coords[[row, 1]], l = plt_coords[[row, 2]])
#   if (row == 1) {
#     areas <- a
#   } else {
#     areas <- c(areas, a)
#   }
# }


plt_dim <- length(venn_plts) %>% sqrt() %>% ceiling()
venn_plts <- wrap_plots(
  venn_plts,
  ncol = plt_dim,
  widths = 2,
  heights = 1
)

ggsave(
  filename = "cell_type_shared_genes_venn.png",
  plot = venn_plts,
  path = args$out_dir,
  width = plt_dim * 4,
  height = plt_dim * 2,
  units = "in",
  dpi = 300
)

# ------------------------------------------------------------------------------

message("\nSESSION INFO\n")
print(sessionInfo())

