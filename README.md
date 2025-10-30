This repository contains code for the manuscript:

Integration of human and mouse single-cell transcriptomes of the developing cerebellum identifies potential cell-of-origin for Group 3 and 4 medulloblastoma
Cheong, I., Lau L., and S. Pai. 
*bioRXiv preprint link TBA*

# All figures in this manuscript
Code to generate all figures in this manuscript in 

# Preparation of input datasets
See `scrnaseq_Leo/integrations/pipeline_instructions.txt` for steps.
Code to standardize processing of each input dataset is in: `scrnaseq_Leo/dataset_processing`

# Human-mouse integration
* The pipeline to perform integrations is in `scrnaseq_Leo/integrations`. See `pipeline_instructions.txt`.
* UMAP of resulting integration: `manuscript_output/figures/results/integ_full.qmd`
* Visual comparison of different integration methods: `manuscript_output/figures/results/integ_methods.qmd`
* Assignment of labels to cell clusters: `software/utilities/cell_labeling.R`
* Subsetting of rhombic lip cells and non-neuronal reference cells: `scrnaseq_Leo/notebook/notes/_20240524.qmd`
* Plotting UMAP, cluster markers of integrated rhombic lip cells: `manuscript_output/figures/results/rl_gene_expr.qmd` and `integ_rl.qmd`

# Subclustering of UBC and plots
* Subsetting for UBC: `lab_notebook/notes/_20240825.qmd`
* Finding and plotting best clustering resolution: `manuscript_output/figures/results/ubc_cluster_res.qmd`
* Plotting UMAP and proportion tests: `manuscript_output/figures/results/ubc_cluster_res.qmd`
* scProportionTest calculate: `scrnaseq_leo/notebook/notes/_20240718.qmd`
* Propeller calculate: `scrnaseq_leo/notebook/notes/_20240606.qmd`
* DEG analysis: `scrnaseq_leo/notebook/notes/_20240715.qmd`
* Pathway analysis: `lab_notebook/notes/_20241029.qmd`
* PySCENIC: `lab_notebook/notes/_20240704.qmd`
* Plots of DEG in UBC clusters and top regulons: `manuscript_output/figures/results/ubc_genes_regulons.qmd`

# Integration with tumours
* QC & testing multiple integration methods for tumours, inc. FastMNN: `lab_notebook/notes/_20230510.qmd`
* Plotting tumour UMAP, Group 4 MB tumour volcano plot, intersection of TFs in Group 4 MB and UBC 1, and top regulons: `manuscript_output/figures/results/mb_similarity.qmd`
* DEG analysis: `lab_notebook/notes/_20240827.qmd`
* PySCENIC: `lab_notebook/notes/_20240928.qmd`
* Projection of developmental cell states: `manuscript_output/figures/results/annot_cluster_single_IanUBCs.R`