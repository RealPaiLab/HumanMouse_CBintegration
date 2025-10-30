# ==============================================================================
# Supplementary figure of manuscript showing the distribution of age and
# datasets in the UBCs. This script should be run using the scrnaseq_env conda
# environment.
# ==============================================================================

library(tidyverse)
library(patchwork)
library(Seurat)

source("./utils.R")
source("../../software/utilities/plotting.R")

out_dir <- file.path(
  "/.mounts/labs/pailab/private/projects/HumanMouseUBC/figures/SuppFig_ubc_age_dataset",
  format(Sys.Date(), "%Y%m%d")
)
if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

# colour palette
my_pals <- get_custom_pals()

# ------------------------------------------------------------------------------
# distribution of ages

# load Seurat object for integrated dataset
srat_qs <- get_srat_paths()
ubc_srat <- load_srat(srat_qs["ubc"]) %>%
  pluck("ubc")

# clean up metadata
ubc_srat$human_age <- fct_drop(ubc_srat$human_age)
ubc_srat$mouse_age <- fct_drop(ubc_srat$mouse_age)
ubc_srat$combined_age <- fct_drop(ubc_srat$combined_age)
ubc_srat$dataset_name <- str_remove(ubc_srat$dataset_name, "full_cerebellum_") %>%
  str_replace(
    pattern = "_(.*)",
    replacement = " (\\1)"
  )
ubc_srat$subclusters <- str_replace(
  # convert "UBC_0" to "iUBC0" etc.
  string = ubc_srat$subclusters,
  pattern = "UBC_([:digit:]{1})",
  replacement = "iUBC\\1"
)

# get ages for which there are <100 human cells
human_ages_to_filter <- table(ubc_srat$human_age) %>%
  as.data.frame() %>%
  filter(Freq < 100) %>%
  pull(1) %>%
  as.character()

# get ages for which there are <50 mouse cells
mouse_ages_to_filter <- table(ubc_srat$mouse_age) %>%
  as.data.frame() %>%
  filter(Freq < 50) %>%
  pull(1) %>%
  as.character()

# get proportion of cells per age
combined_age <- table(ubc_srat$subclusters, ubc_srat$combined_age, ubc_srat$species) %>%
  as.data.frame() %>%
  rename(subcluster = Var1, age = Var2, species = Var3, num_cells = Freq) %>%
  filter(!age %in% c(human_ages_to_filter, mouse_ages_to_filter) & num_cells != 0)

# plot distribution of ages
age_distribution <- ggplot(
  data = combined_age,
  mapping = aes(
    x = age,
    y = num_cells,
    fill = subcluster,
    group = subcluster
  )
) + 
  geom_area(position = "fill") + 
  facet_wrap(
    vars(species),
    ncol = 2,
    scales = "free_x",
    labeller = labeller(species = str_to_title)
  ) + 
  labs(
    caption = sprintf(
      "Filtered human cells from %s PCW, mouse cells from %s",
      paste0(str_remove(human_ages_to_filter, pattern = " PCW"), collapse = "/"),
      paste0(mouse_ages_to_filter, collapse = "/")
    ),
    y = "proportion",
    fill = "UBC\nsubcluster"
  ) + 
  scale_fill_manual(values = my_pals$ubc_integ_clust2) + 
  scale_x_discrete(expand = expansion(mult = 0)) + 
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) + 
  theme_classic2() + 
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1),
    panel.spacing = unit(0.05, "npc")
  )

ggsave(
  filename = "age_distribution.pdf",
  plot = age_distribution,
  path = out_dir,
  width = 9,
  height = 6,
  units = "in"
)

# ------------------------------------------------------------------------------
# distribution of datasets

# UMAP coloured by dataset
dataset_umap <- DimPlot(
  ubc_srat,
  reduction = "umap",
  group.by = "dataset_name",
  label = FALSE
) + 
  labs(title = "Human/mouse UBCs", x = "UMAP 1", y = "UMAP 2") + 
  theme_classic2() + 
  theme(legend.position = "none")
ggsave(
  filename = "dataset_umap.pdf",
  plot = dataset_umap,
  path = out_dir,
  width = 4,
  height = 4,
  units = "in"
)

# proportion of cells from each dataset in each cluster
dataset_prop <- cluster_barplot(
  object = ubc_srat,
  split.by = "dataset_name",
  group.by = "subclusters",
  position = "fill"
) + 
  labs(x = "UBC subcluster", fill = "dataset") + 
  theme_classic2() + 
  theme(axis.text.x = element_text(angle = 30, hjust = 1))
ggsave(
  filename = "dataset_prop_bar.pdf",
  plot = dataset_prop,
  path = out_dir,
  width = 5,
  height = 4,
  units = "in",
  dpi = 600
)

# ------------------------------------------------------------------------------
# final figure

layout <- c(
  # age plot
  area(t = 1, b = 3, l = 1, r = 5),
  # dataset UMAP
  area(t = 4, b = 5, l = 1, r = 2),
  # dataset bar plot
  area(t = 4, b = 5, l = 3, r = 5)
)

.plt <- wrap_plots(
  free(age_distribution) + labs(caption = NULL),
  dataset_umap,
  free(dataset_prop)
) + 
  plot_layout(design = layout)

walk(
  .x = c("png", "pdf"),
  .f = \(device) {
    ggsave(
      filename = paste0("combined.", device),
      plot = .plt,
      path = out_dir,
      width = 7.5,
      height = 6.5,
      units = "in",
      dpi = 600
    )
  }
)
