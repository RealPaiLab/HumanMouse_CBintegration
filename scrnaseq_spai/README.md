# Human UBC clustering and analysis

* Whole cerebellum: [CCA integration of Aldinger and Sepp](https://github.com/RealPaiLab/MB_scRNAseq/blob/spai/scrnaseq_Leo/integrations/HPC_humanCBOnly.sh) - note: this code is in `scrnaseq_Leo` of `spai` branch
* Whole cerebellum: [Assign cluster identity](https://github.com/RealPaiLab/MB_scRNAseq/blob/spai/scrnaseq_spai/AssignClusterIdentity.R) and subset for rhombic lip
* Start with RL and [subset for UBC](https://github.com/RealPaiLab/MB_scRNAseq/blob/spai/scrnaseq_spai/profileRL_reclusterUBC.R)
* Subset for EOMES+ UBC and cluster: [Code with Quang](https://github.com/RealPaiLab/MB_scRNAseq/blob/qtrinh/scrnaseq_qtrinh/scripts/UBC.main.R)
* Human UBC clusters: [Differential expression analysis](https://github.com/RealPaiLab/MB_scRNAseq/blob/spai/scrnaseq_spai/humanUBCclusters/FromIan/01A_diff_expr.R) and heatmap, volcano plots
* **Ian's code:** [Characterizing human-specific UBC clusters](https://github.com/RealPaiLab/MB_scRNAseq/tree/analysis_ian/HumanMouseUBC/integrated_human_ubc)

  **NOTE:** CytoTRACE, pySCENIC Scenic+ runs on the cluster. You will need special conda environments for this.
  Speak to Shraddha about getting these. 
  

# Human/mouse integration

- CCA/RPCA/Harmony integration
  - [job submission scripts](https://github.com/RealPaiLab/MB_scRNAseq/tree/analysis_ian/scrnaseq_Leo/integrations)
  - [pipeline scripts](https://github.com/RealPaiLab/MB_scRNAseq/tree/analysis_ian/scrnaseq_Leo/pipeline_scripts) (if more detail is needed)
- clustering resolutions: [clustree](https://github.com/RealPaiLab/MB_scRNAseq/blob/analysis_ian/scrnaseq_Leo/notebook/notes/_20240527.qmd)
- Testing for proportions
  - [scProportionTest](https://github.com/RealPaiLab/MB_scRNAseq/blob/analysis_ian/scrnaseq_Leo/notebook/notes/_20240718.qmd)
  - [propeller](https://github.com/RealPaiLab/MB_scRNAseq/blob/analysis_ian/lab_notebook/notes/_20240606.qmd)
- Differential gene expression: [see Leo's notes](https://github.com/RealPaiLab/MB_scRNAseq/blob/analysis_ian/scrnaseq_Leo/notebook/notes/_20240715.qmd)
- Pathway enrichment analysis: [g:Profiler](https://github.com/RealPaiLab/MB_scRNAseq/blob/analysis_ian/lab_notebook/notes/_20241029.qmd)
- Notes on [generation of custom GMT files](https://github.com/RealPaiLab/MB_scRNAseq/blob/analysis_ian/lab_notebook/notes/_paths_to_datasets.qmd). (See "Bader lab gene sets" section.)
As described in the notes, the actual GMT files are located in `/.mounts/labs/pailab/src/BaderLab_EM_Genesets/August_08_2023/Human/symbol/`. I filtered the GMT file we got from the EnrichmentMap website for pathways with 10-250 genes (inclusive). This is saved in `Human_GOBP_AllPathways_no_GO_iea_August_08_2023_symbol_min10_max250.gmt`. I then uploaded this new GMT file to g:Profiler, which returns an ID - this ID is what you need to use as the "organism" parameter when running `gprofiler2::gost()` in R. There's a table with the IDs in the notes I sent.
- Pseudotime: [Monocle 3](https://github.com/RealPaiLab/MB_scRNAseq/blob/analysis_ian/software/integrated/run_monocle.R)
- GRN inference: [pySCENIC](https://github.com/RealPaiLab/MB_scRNAseq/tree/analysis_ian/software/pyscenic)
  - [converting Seurat to CSV pySCENIC input](https://github.com/RealPaiLab/MB_scRNAseq/blob/analysis_ian/software/pyscenic/seurat_to_scenic.R)
  - pySCENIC: Jupyter notebook in `results/pyscenic/20240704/` (in isilon only, file was too big to commit to GitHub)

# Vladoiu et al. MB tumours

- Integration/batch correction: [FastMNN](https://github.com/RealPaiLab/MB_scRNAseq/blob/analysis_ian/software/tumour/Vladoiu/test_integ_methods.R)
- GRN inference: [pySCENIC](https://github.com/RealPaiLab/MB_scRNAseq/tree/analysis_ian/software/pyscenic)
  - [converting Seurat to CSV pySCENIC input](https://github.com/RealPaiLab/MB_scRNAseq/blob/analysis_ian/software/pyscenic/seurat_to_scenic.R)
  - pySCENIC: Jupyter notebook in `results/pyscenic/20231018/` (in isilon only, file was too big to commit to GitHub)
 
# Ian's Figures

- [Thesis figures](https://github.com/RealPaiLab/MB_scRNAseq/tree/analysis_ian/thesis/figures)
- [Figures for various presentations/posters](https://github.com/RealPaiLab/MB_scRNAseq/tree/analysis_ian/figures)
- [Human/mouse UBC manuscript figures (work in progress)](https://github.com/RealPaiLab/MB_scRNAseq/tree/analysis_ian/HumanMouseUBC/figures)


