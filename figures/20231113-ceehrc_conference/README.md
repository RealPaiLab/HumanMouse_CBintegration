## Regulon expression

Regulon expression CSV files were generated from:

-   `/CBL_scRNAseq/results/pyscenic/20230725/aldinger_RL.h5ad`
-   `/CBL_scRNAseq/results/pyscenic/20231018/vladoiu_mb.h5ad`

The data were extracted from h5ad in Python, similar to below:

```python
import anndata as ad

adata = ad.read_h5ad('<input h5ad file>')
adata.obs.to_csv('<output csv file>')
```

This is essentially the same as a Seurat object's metadata (`srat[[]]`).


## DepMap plot

Code was copied from [Shraddha's GitHub](https://github.com/RealPaiLab/FetalHindbrain_Epigenetics/blob/master/CRISPR_target_prioritization/prioritizeCandidates.R) and modified.
