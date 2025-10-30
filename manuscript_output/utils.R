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
    merge = "/.mounts/labs/pailab/private/icheong/CBL_scRNAseq/results/integrated/20241223/merged_seurat.qs"
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
      if (endsWith(x = pth, suffix = "qs")) {
        qs::qread(pth)
      } else if (endsWith(x = pth, suffix = "rds")) {
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


# ------------------------------------------------------------------------------
# the functions below were provided by Xinghan (`get_neurodev_g34_genes_func.R`)

# (from Xinghan)
#' get neurodevelopment genes
#' @param path (charaters) The path to the file containing the gene list
#' @return (character) A vector of neurodevelopment genes
get_neurodevGenes <- function() {
  additional_genes = c("OLIG3") # mentioned by Kim Aldinger
  
  paths <- c(avc = paste0("/.mounts/labs/pailab/src/gene_list/brain-development",
                        "/aldinger2021_vladoiu2019_carter2018.csv"),
             leto = paste0("/.mounts/labs/pailab/src/gene_list/brain-development",
                  "/LetoEtAl.txt"),
             bhaduri = paste0("/.mounts/labs/pailab/src/gene_list/brain-development",
                            "/BhaduriKriegstein2021_nature_neocortex.txt"),
             sepp = paste0("/.mounts/labs/pailab/src/gene_list/brain-development",
                         "/Sepp2023_CB.txt"),
             ian_leo = paste0("/.mounts/labs/pailab/src/gene_list/brain-development",
                            "/cell_gene_mapping.csv")
             )
  
  # aldinger2021_vladoiu2019_carter2018.csv #
  avc <- read.csv(paths["avc"], stringsAsFactors = F, header = T)
  # select only stem/neuro progenitors
  avc_genes <- unique(avc[grepl("RL|UBCs|GCPs|stem", avc$region), "human"])
  
  # LetoEtAl.txt #
  leto <- read.table(paths["leto"], stringsAsFactors = F, header = F, sep = "\t")
  leto_genes <- unique(leto$V2)
  
  # BhaduriKriegstein2021_nature_neocortex.txt #
  bhaduri <- read.table(paths["bhaduri"], 
                        stringsAsFactors = F, header = F, sep = "\t"
                        )
  bhaduri_genes <- bhaduri[! grepl("OPC|astrocytes", bhaduri$V2), "V1"]
  
  # Sepp2023 #
  sepp <- read.table(paths["sepp"], stringsAsFactors = F, header = F, sep = "\t")
  sepp_genes <- unlist(stringr::str_split(sepp$V2, ", "))
  sepp_genes <- stringr::str_trim(sepp_genes)
  sepp_genes <- unique(sepp_genes)
  
  # Ian & Leo #
  ianLeo <- read.csv(paths["ian_leo"], stringsAsFactors = F, header = T)
  ianLeo <- ianLeo[!ianLeo$reference %in% c("DEA", ""),]
  ianLeo <- ianLeo[!grepl("(oligodendrocytes|Microglia|Endothelial|Astrocytes|Purkinje)", 
                         ianLeo$lineage),]
  ianLeo_genes <- unique(ianLeo$gene)
  
  ### Combine ###
  res <- unique(c(avc_genes, leto_genes, bhaduri_genes, sepp_genes, ianLeo_genes,
                  additional_genes)
                )
  message(sprintf("Obtained %i neuro development genes from %s", 
                  length(res),
                  paste(c("", paths), collapse = "\n- ")
  )
  )
  
  return(res)
}


# (from Xinghan)
#' get G3/4 MB genes
#' @param db (character) The directory to MB_gene list
#' @return (character) A vector of Grp3/4 MB genes
get_g34genes <- function(db = paste0("/.mounts/labs/pailab/src/gene_list/MB_gene",
                                   "/MBgene_database_20240709171001.csv")
                         ) {
  df <- read.csv(db, stringsAsFactors = F, header = T, row.names = 1)
  sub <- df[,grepl("Hendrikse2022_G3_G4_MB_genes|Northcott2017_G34_genes", 
                   colnames(df)
                   )
            ]

  g34_genes <- unique(rownames(sub[rowSums(sub) != 0,]))
  
  return(g34_genes)
}
