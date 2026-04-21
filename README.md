This repository contains code for the manuscript: 
**Integration of human and mouse single-cell transcriptomes of the developing cerebellum identifies potential cell-of-origin for Group 3 and 4 medulloblastoma**

Ian Cheong, Leo Lau, Shraddha Pai

[Read the paper](https://link.springer.com/article/10.1186/s12885-026-16006-1)

## Preparation of input datasets
See `scrnaseq_Leo/integrations/pipeline_instructions.txt` for steps.
Code to standardize processing of each input dataset is in: `scrnaseq_Leo/dataset_processing`

## Human-only analysis
The pipeline to perform integrations is in `scrnaseq_Leo/integrations`. See `HPC_humanCBOnly.sh`. The rest of the characterization code is in `manuscript_output/humanUBCclustering`
* Clustering UBC: `HumanOnly_FindBestUBCCluster.R`
* Differential expression analysis: `01A_diff_expr.R`
* Pathway enrichement: `01B_pathway_enrichment.R`
* Cytoscape viz: `01C_cytoscape.R`
* Cytotrace: `02A_cytotrace.R`
* UBC clustering for Aldinger & Sepp PAX6+ EOMES+ cells separately, and consensus clusters: `humanSingleDatasetUBCclustering.R`

## Mouse-only analysis:
In the `manuscript_output/figures/results` directory
* Mouse UBC clusters: `mouseUBC_clusters.R`
* Identify clusters in Vladoiu and Sepp datasets separately, and find consensus: `mouseSingleDatasetUBCclusters.R`

## Human-mouse integration
* The pipeline to perform integrations is in `scrnaseq_Leo/integrations`. See `pipeline_instructions.txt`.
* UMAP of resulting integration: `manuscript_output/figures/results/integ_full.qmd`
* Visual comparison of different integration methods: `manuscript_output/figures/results/integ_methods.qmd`
* Assignment of labels to cell clusters: `software/utilities/cell_labeling.R`
* Subsetting of rhombic lip cells and non-neuronal reference cells: `scrnaseq_Leo/notebook/notes/_20240524.qmd`
* Plotting UMAP, cluster markers of integrated rhombic lip cells: `manuscript_output/figures/results/rl_gene_expr.qmd` and `integ_rl.qmd`
* Proportion tests for rhombic lip: `lab_notebook/notes/_20240822.qmd` and `lab_notebook/notes/_20240823.qmd`

## Subclustering of UBC and plots
* Subsetting for UBC: `lab_notebook/notes/_20240825.qmd`
* Finding and plotting best clustering resolution: `manuscript_output/figures/results/ubc_cluster_res.qmd`
* Plotting UMAP and proportion tests: `manuscript_output/figures/results/ubc_cluster_res.qmd`
* scProportionTest calculate: `scrnaseq_leo/notebook/notes/_20240718.qmd`
* Propeller calculate: `lab_notebook/notes/_20240606.qmd`
* DEG analysis: `scrnaseq_leo/notebook/notes/_20240715.qmd`
* Pathway analysis: `lab_notebook/notes/_20241029.qmd`
* PySCENIC: `lab_notebook/notes/_20240704.qmd`
* Plots of DEG in UBC clusters and top regulons: `manuscript_output/figures/results/ubc_genes_regulons.qmd`
* Pseudotime calculation: `lab_notebook/notes/_20240611.qmd`
* Pseudotime plots: `manuscript_output/figures/results/pseudotime.qmd`
* Plot cluster distribution by sex: `manuscript_output/figures/results/plot_ubc_sexdistr.qmd`
* Milo analysis: `manuscript_output/figures/results/miloProp.R`

## Integration with tumours
* QC & testing multiple integration methods for tumours, inc. FastMNN: `lab_notebook/notes/_20230510.qmd`
* Plotting tumour UMAP, Group 4 MB tumour volcano plot, intersection of TFs in Group 4 MB and UBC 1, and top regulons: `manuscript_output/figures/results/mb_similarity.qmd`
* DEG analysis: `lab_notebook/notes/_20240827.qmd`
* PySCENIC: `lab_notebook/notes/_20240928.qmd`
* Projection of developmental cell states: `manuscript_output/figures/results/annot_cluster_single_IanUBCs.R`
