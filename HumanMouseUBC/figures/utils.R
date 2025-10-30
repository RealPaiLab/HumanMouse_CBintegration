# ==============================================================================
# Convenience functions for thesis figures.
# ==============================================================================

#' Get paths to Seurat objects.
#'
#' @return Named character vector containing with the full path to the Seurat
#'   objects for: Vladoiu mouse dataset, Aldinger human dataset, Sepp human
#'   dataset, Sepp mouse dataset, integrated full cerebellum, integrated RL
#'   lineage, integrated UBCs.
#'
get_srat_paths <- function() {
  # Seurat file paths for each individual dataset
  srat_qs <- read_csv(
    "/.mounts/labs/pailab/private/icheong/CBL_scRNAseq/scrnaseq_Leo/integrations/dataset_list.csv"
  ) %>%
    select(-dataset_location_rds) %>%
    filter(
      str_detect(dataset_name, "_full_cerebellum_"),
      !str_detect(dataset_name, "Luo|Vladoiu")
    ) %>%
    mutate(
      dataset_name = tolower(str_remove(dataset_name, "full_cerebellum_"))
    ) %>%
    tibble::deframe()

  # Seurat file paths for the integrated datasets
  srat_qs <- c(
    srat_qs,
    # Vladoiu dataset (with UMAP)
    vladoiu_mouse = "/.mounts/labs/pailab/private/icheong/CBL_scRNAseq/results/mouse/Vladoiu/20240810/seurat_with_umap.qs",
    # full cerebellum - CCA
    full = "/.mounts/labs/pailab/private/llau/results/integrated/20240516/cca/20240516_cca_integ.qs",
    # RL lineage (+ controls) - CCA
    rl = "/.mounts/labs/pailab/private/llau/results/integrated/20240524/25_pc_without_luo/25_pc_rl.qs",
    # RL lineage (+ controls) - RPCA
    rl_rpca = "/.mounts/labs/pailab/private/llau/results/integrated/20240613/rpca/combined_clusters_rl.qs",
    # RL lineage (+ controls) - harmony
    rl_harmony = "/.mounts/labs/pailab/private/llau/results/integrated/20240604/harmony_rl/rl.qs",
    # UBC subclusters
    ubc = "/.mounts/labs/pailab/private/icheong/CBL_scRNAseq/results/integrated/20240825/ubc_subset.qs",
    # medulloblastoma dataset from Vladoiu
    mb = "/.mounts/labs/pailab/private/icheong/CBL_scRNAseq/results/tumour/Vladoiu/20230510/mb_mnn.rds",
    # full cerebellum - merge only (no integration)
    merge = "/.mounts/labs/pailab/private/icheong/CBL_scRNAseq/results/integrated/20241223/merged_seurat.qs",
    # human UBC clusters only (from Quang)
    quang_ubc = "/.mounts/labs/pailab/private/projects/HumanMouseUBC/UBCclusters/fromQuang/UBC.Harmony.RDS"
  )

  return(srat_qs)
}


#' Read in Seurat objects from a `qs` file.
#'
#' @param path Vector of path(s) to the Seurat objects to read in.
#'
#' @return List of Seurat objects.
#'
load_srat <- function(
  path = get_srat_paths()
) {
  srat <- map2(
    .x = path,
    .y = names(path),
    .f = \(pth, nme) {
      message(sprintf("***Reading in %s from %s***", nme, pth))
      if (endsWith(tolower(pth), "qs")) {
        qs::qread(pth)
      } else if (endsWith(tolower(pth), "rds")) {
        readRDS(pth)
      }
    }
  )
}


#' Get custom colour palettes.
#'
#' @return Named list of colour palettes: `species` (for human/mouse),
#'   `rl_integ_clust` (for CCA-integrated RL lineage clusters), `rl_integ_annot`
#'   (for CCA-integrated broad RL lineage cell types)
#'
get_custom_pals <- function() {
  custom_pals <- list(
    # colour palette for species
    species = c("grey", scales::pal_hue()(1)),
    # colour palette for CCA-integrated RL lineage clusters
    rl_integ_clust = pals::cols25(18),
    # colour palette for CCA-integrated RL lineage broad cell types
    rl_integ_annot = setNames(
      # remove yellow, too similar to the UBC yellow
      pals::trubetskoy(8)[-3] %>% unname(),
      nm = c("endothelial", "GC", "GCP", "microglia", "oligodendrocyte/OPC", "RL", "UBC")
    ),
    # colour palette for CCA-integrated UBC clusters
    ubc_integ_clust = setNames(
      pals::brewer.set2(6),
      nm = paste0("UBC_", c(0:5))
    ),
    # colour palette for CCA-integrated UBC clusters (same as above, but with
    # different names)
    ubc_integ_clust2 = setNames(
      pals::brewer.set2(6),
      nm = paste0("iUBC", c(0:5))
    ),
    mb_subtype = setNames(
      pals::okabe(8)[c(8, 3, 4)],
      nm = c("SHH", "G3", "G4")
    )
  )
  return(custom_pals)
}


#' Custom `ggplot2` theme based on `theme_classic()`.
#'
#' @param base_size Base font size in points
#' @param base_family Base font family
#' @param base_line_size Base size for line elements
#' @param base_rect_size Bsae size for rect elements
#'
#' @return A `ggplot2` theme.
#'
theme_classic2 <- function(
  base_size = 11,
  base_family = "",
  base_line_size = base_size / 22,
  base_rect_size = base_size / 22
) {
  theme_classic(
    base_size = base_size,
    base_family = base_family,
    base_line_size = base_line_size,
    base_rect_size = base_rect_size
  ) %+replace% 
    theme(
      axis.text = element_text(colour = "black", size = rel(0.8)),
      axis.ticks = element_line(colour = "black"),
      legend.key.size = unit(0.8, "lines"),
      legend.box.spacing = unit(0.25 * base_size, "points"),
      plot.title = element_text(size = rel(1.2), hjust = 0.5),
      plot.tag = element_text(face = "bold", hjust = 0.5, vjust = 0.5),
      plot.tag.location = "plot",
      strip.background = element_blank(),
      strip.text = element_text(size = rel(1)),
      complete = TRUE
    )
}


add_plot_tags <- function(
  plt_list,
  tags
) {
  plt_list <- map2(
    .x = plt_list,
    .y = tags,
    .f = \(plt, tag) {
      plt <- plt + 
        labs(tag = tag)
    }
  )
}

