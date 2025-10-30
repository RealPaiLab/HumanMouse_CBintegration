# List of scripts

These scripts were used to analyze the integrated human and mouse scRNA-seq data.

## `diff_exp_genes.R`

Takes integrated Aldinger/Vladoiu Seurat object - cluster 8 contains both human and mouse UBCs while clusters 17/18 contain mostly human UBCs. This script identifies differentially expressed genes between human UBCs from cluster 8 and human UBCs from clusters 17/18. In other words, what's the difference between the human UBCs which cluster with other mouse UBCs and those that don't?

## `integrate_liamHs_vladoiuMm.R` and `visualize_liamHs_vladoiuMm.R`

Integrate the human RL (Seurat object from Liam) with the entire mouse cerebellum (data from Vladoiu). Performs clustering and dimensionality reduction and saves the resulting Seurat object. Creates various plots for visualization.

## `integrate_liamHs_vladoiuMmRL.R` and `visualize_liamHs_vladoiuMmRL.R`

Integrate the human RL (Seurat object from Liam) with the mouse RL (data from Vladoiu). Performs clustering and dimensionality reduction and saves the resulting Seurat object. Creates various plots for visualization.
