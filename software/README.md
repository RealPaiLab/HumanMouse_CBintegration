
# Analyses

## Human/mouse integration

- CCA/RPCA/Harmony integration
  - [job submission scripts](https://github.com/RealPaiLab/MB_scRNAseq/tree/analysis_ian/scrnaseq_Leo/integrations)
  - [pipeline scripts](https://github.com/RealPaiLab/MB_scRNAseq/tree/analysis_ian/scrnaseq_Leo/pipeline_scripts) (if more detail is needed)
- clustering resolutions: [clustree](https://github.com/RealPaiLab/MB_scRNAseq/blob/analysis_ian/scrnaseq_Leo/notebook/notes/_20240527.qmd)
- Testing for proportions
  - [scProportionTest](https://github.com/RealPaiLab/MB_scRNAseq/blob/analysis_ian/scrnaseq_Leo/notebook/notes/_20240718.qmd)
  - [propeller](https://github.com/RealPaiLab/MB_scRNAseq/blob/analysis_ian/lab_notebook/notes/_20240606.qmd)
- Differential gene expression: [see Leo's notes](https://github.com/RealPaiLab/MB_scRNAseq/blob/analysis_ian/scrnaseq_Leo/notebook/notes/_20240715.qmd)
- Pathway enrichment analysis: [g:Profiler](https://github.com/RealPaiLab/MB_scRNAseq/blob/analysis_ian/lab_notebook/notes/_20241029.qmd)
- Pseudotime: [Monocle 3](https://github.com/RealPaiLab/MB_scRNAseq/blob/analysis_ian/software/integrated/run_monocle.R)
- GRN inference: [pySCENIC](https://github.com/RealPaiLab/MB_scRNAseq/tree/analysis_ian/software/pyscenic)
  - [converting Seurat to CSV pySCENIC input](https://github.com/RealPaiLab/MB_scRNAseq/blob/analysis_ian/software/pyscenic/seurat_to_scenic.R)
  - pySCENIC: Jupyter notebook in `results/pyscenic/20240704/` (in isilon only, file was too big to commit to GitHub)

## Vladoiu et al. MB tumours
- SingleR projection from development: `manuscript_figures/results/annot_cluster_single_IanUBC.R`
- Integration/batch correction: [FastMNN](https://github.com/RealPaiLab/MB_scRNAseq/blob/analysis_ian/software/tumour/Vladoiu/test_integ_methods.R)
- GRN inference: [pySCENIC](https://github.com/RealPaiLab/MB_scRNAseq/tree/analysis_ian/software/pyscenic)
  - [converting Seurat to CSV pySCENIC input](https://github.com/RealPaiLab/MB_scRNAseq/blob/analysis_ian/software/pyscenic/seurat_to_scenic.R)
  - pySCENIC: Jupyter notebook in `results/pyscenic/20231018/` (in isilon only, file was too big to commit to GitHub)
