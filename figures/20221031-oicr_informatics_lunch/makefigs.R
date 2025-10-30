# ==============================================================================
# These figures were generated for the OICR Informatics Lunch and Learn on
# Monday, October 31, 2022.
# ==============================================================================

# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(tidyverse)
library(Seurat)

# import data
srat_names <- c("aldinger", "vladoiu", "cca", "rpca", "harmony")
srat_rds <- c(
  # Aldinger human RL from Liam
  "/isilon/CBL_scRNAseq/data/human/Aldinger/glutamatergic_dev_Liam.RDS",
  # Vladoiu mouse RL
  "/CBL_scRNAseq/results/mouse/Vladoiu/merged_seurat_RLonly.rds",
  # CCA
  "/CBL_scRNAseq/results/integrated/vladoiu_liam_RL.rds",
  # RPCA
  "/CBL_scRNAseq/results/integrated/20221002/kanchors80/vladoiu_liam_RL_rpca.rds",
  # harmony
  "/CBL_scRNAseq/results/integrated/20221003/vladoiu_liam_RL_harmony.rds"
)

all_srat <- purrr::map(
  .x = srat_rds,
  .f = readRDS
)
names(all_srat) <- srat_names
integ_methods <- c("cca", "rpca", "harmony")

print(all_srat)

# import functions
root_dir <- "CBL_scRNAseq"
source(file.path("", root_dir, "software/mouse/Vladoiu/add_annotations.R"))
source(file.path("", root_dir, "software/utilities/cell_labelling.R"))
source(file.path("", root_dir, "software/utilities/score_integration.R"))

# ------------------------------------------------------------------------------
# separate UMAPs of human and mouse

# human aldinger
.plt <- DimPlot(
  all_srat$aldinger,
  reduction = "umap",
  group.by = "new_cell_type",
  label = TRUE,
  repel = TRUE
) + 
  NoLegend() + 
  labs(title = NULL)

ggsave(
  "umap_aldinger.png",
  plot = .plt,
  width = 5,
  height = 5,
  units = "in",
  dpi = 600
)

# mouse vladoiu
.plt <- DimPlot(
  all_srat$vladoiu,
  reduction = "umap",
  group.by = "seurat_clusters",
  label = TRUE,
  repel = TRUE
) + 
  NoLegend() + 
  labs(title = NULL)

ggsave(
  "umap_vladoiu.png",
  plot = .plt,
  width = 5,
  height = 5,
  units = "in",
  dpi = 600
)

# CCA integration
.plt <- DimPlot(
  all_srat$cca,
  reduction = "umap",
  group.by = "seurat_clusters",
  label = TRUE,
  repel = TRUE
) + 
  NoLegend() + 
  labs(title = NULL)

ggsave(
  "umap_cca_clusters.png",
  plot = .plt,
  width = 5.5,
  height = 5,
  units = "in",
  dpi = 600
)

# ------------------------------------------------------------------------------
# UMAPs of integrated data

# colour palette
species_cols <- c(hcl(h = 15, c = 100, l = 65), "grey")

for (srat in integ_methods) {
  .plt <- DimPlot(
    all_srat[[srat]],
    label = FALSE,
    group.by = "species",
    pt.size = 0.01,
    cells.highlight = rownames(all_srat[[srat]]@meta.data)[all_srat[[srat]]$species == "human"],
    sizes.highlight = 0.01
  ) + 
    NoLegend() + 
    scale_color_manual(labels = c("mouse", "human"), values = rev(species_cols)) + 
    labs(title = NULL)
  
  ggsave(
    paste0("umap_", srat, "_species.png"),
    plot = .plt,
    width = 5.5,
    height = 5,
    units = "in",
    dpi = 600
  )
  
  # add mouse cell types (`add_annotation.R`)
  all_srat[[srat]] <- label_vladoiu_cells(all_srat[[srat]]) %>% 
    # pool cell types (`pool_cell_types.R`)
    pool_cell_types(col_name = "common_cell_type")
  
  # change "other/missing" to NA
  all_srat[[srat]]@meta.data <- mutate(
    all_srat[[srat]]@meta.data,
    common_cell_type = case_when(
      common_cell_type != "other/missing" ~ common_cell_type
    )
  )
  
  .plt <- DimPlot(
    all_srat[[srat]],
    group.by = "common_cell_type",
    label = TRUE,
    repel = TRUE
  ) + 
    labs(title = NULL)
  
  ggsave(
    paste0("umap_", srat, "_celltype.png"),
    plot = .plt,
    width = 6,
    height = 5,
    units = "in",
    dpi = 600
  )
}

# ------------------------------------------------------------------------------
# EOMES and LMX1A expression in CCA

for (gene in c("EOMES", "LMX1A")) {
  .plt <- FeaturePlot(
    all_srat$cca,
    features = gene,
    min.cutoff = "q10",
    max.cutoff = "q90",
    order = TRUE
  )
  
  ggsave(
    paste0(gene, ".png"),
    plot = .plt,
    width = 5,
    height = 4,
    units = "in",
    dpi = 600
  )
}

# ------------------------------------------------------------------------------
# bar chart of number of human/mouse cells per cluster

num_cells <- table(all_srat$cca$species, all_srat$cca$seurat_clusters) %>% 
  as.data.frame(.)
colnames(num_cells) <- c("species", "cluster", "freq")

# plot absolute number of human/mouse cells in each cluster
.plt <- ggplot(num_cells, aes(x = cluster, y = freq, fill = fct_rev(species))) + 
  geom_col() + 
  labs(y = "number") + 
  theme_light() + 
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) + 
  scale_fill_manual(values = rev(species_cols), name = "species") + 
  theme_classic() + 
  theme(legend.position = c(0.9, 0.9), 
        axis.text = element_text(colour = "black"), 
        axis.text.x = element_text(angle = 45, hjust = 1), 
        axis.ticks = element_line(colour = "black"))

ggsave(
  "species_bar.png",
  plot = .plt,
  width = 8,
  height = 5,
  units = "in",
  dpi = 600
)

# ------------------------------------------------------------------------------
# plot distance between cluster centroids

# get distance between human and mouse centroids
all_centroid_dist <- purrr::map2_dfr(
  .x = all_srat[names(all_srat) %in% integ_methods],
  .y = names(all_srat)[names(all_srat) %in% integ_methods],
  .f = function(x, y) {
    # add mouse cell types (`add_annotation.R`)
    x <- label_vladoiu_cells(x) %>% 
      # pool cell types (`pool_cell_types.R`)
      pool_cell_types(col_name = "common_cell_type")
    
    # get embeddings and metadata
    if ("harmony" %in% names(x@reductions)) {
      embeddings <- Embeddings(x, reduction = "harmony")
    } else {
      embeddings <- Embeddings(x, reduction = "pca")
    }
    metadata <- x@meta.data
    
    # calculate distance between human and mouse centroids for each cell type
    all_dist <- delta_centroids(
      embeddings = embeddings,
      metadata = metadata,
      by = "common_cell_type",
      filter_out = c("other/missing"),
      integ_method = y
    )
    
    return(all_dist)
  }
)

# fix order of integ_method
all_centroid_dist <- mutate(
  all_centroid_dist,
  integ_method = fct_relevel(integ_method, "cca", "rpca", "harmony")
)

# make plot
.plt <- ggplot(all_centroid_dist, aes(x = common_cell_type, y = delta)) + 
  geom_point(aes(colour = integ_method, shape = integ_method), size = 3) + 
  labs(
    x = "Cell type",
    y = "Distance between human\nand mouse cluster centroids",
    colour = "Integration\nmethod",
    shape = "Integration\nmethod"
  ) + 
  scale_y_continuous(limits = c(0, NA)) + # `NA` uses the current min/max
  scale_colour_brewer(palette = "Set1", labels = c("CCA", "RPCA", "Harmony")) + 
  scale_shape_discrete(labels = c("CCA", "RPCA", "Harmony")) + 
  theme_light() + 
  theme(axis.text = element_text(colour = "black"))

ggsave(
  filename = "delta_centroids.png",
  plot = .plt,
  width = 4,
  height = 3,
  units = "in",
  dpi = 600
)
