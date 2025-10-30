## Setting up conda on the HPC

**!!NOTE. These instructions are for setting up conda on the HPC, not the pailab server.**

First, make sure that conda is available in modulator (see the OICR modulator repo for more information). Then, follow the steps below:

```bash
module load conda/24.1.2
conda init

# change default conda directories to isilon - more space available
cd /.mounts/labs/pailab/private/<username>/
mkdir -p .conda/envs .conda/pkgs
conda config --add envs_dirs /.mounts/labs/pailab/private/<username>/.conda/envs/
conda config --add pkgs_dirs /.mounts/labs/pailab/private/<username>/.conda/pkgs/
```

To make sure everything is set up, check the `~/.condarc` file, it should show something like:

```
envs_dirs:
  - /.mounts/labs/pailab/private/icheong/.conda/envs/
pkgs_dirs:
  - /.mounts/labs/pailab/private/icheong/.conda/pkgs/
```

Conda should now be installed and ready to use whenever you log into the cluster (you don't need to load the conda module at all).

## Using the conda YAM files

To set up:
Conda YAML file is located at `MB_scRNAseq/conda`
    a. Create Conda environment by running: 

        conda env create -f <conda_env>.yaml
    
    b. To run environment use:

        conda activate <conda_env>

    c. To update environment in case of new packages or you want to update an existing package

        conda env update --file <conda_env>.yaml --prune

## Contents

This folder contains conda environment files created using `conda env export`.

### `scrnaseq_env`

The default single-cell RNA-seq analysis environment.

Note that the `seurat-wrappers` package must be installed separately. First, activate the `scrnaseq_env` environment:

```bash
# activate the scrnaseq_env environment
$ conda activate scrnaseq_env
```

Then install SeuratWrappers in R (as per the [instructions](https://github.com/satijalab/seurat-wrappers)):

```r
# specify the commit number for reproducibility (also cuz later versions require Seurat v5)
remotes::install_github("satijalab/seurat-wrappers", ref = "d28512f")
```

### `diff_expr`

This environment was created for running differential expression analysis with `FindMarkers` using tests other than the default Wilcoxon rank-sum test. A separate environment was created since `MAST` version 1.26.0 was required; the default `scrnaseq_env` environment could only install `MAST` version 1.24.0 which gave errors when running `FindMarkers`.

The `scrnaseq_env` environment can still be used for `FindMarkers` if `MAST` is not needed.

### `propeller`

This environment was created to test for differences in cell type proportions using `propeller` in the `speckle` package.

### `monocle3`

This environment was created for running Monocle 3 and CytoTRACE.

Note that the actual `monocle3`, `seurat-wrappers`, and `CytoTRACE` packages must be installed separately. First, activate the `monocle3` environment:

#### `monocle3` installation

```bash
# activate the monocle3 environment
$ conda activate monocle3
```

Then install `monocle3` in R (as per the [instructions](https://cole-trapnell-lab.github.io/monocle3/docs/installation/) on the Monocle 3 website):

```r
# specify the commit number for reproducibility
remotes::install_github("cole-trapnell-lab/monocle3", ref = "98402ed", upgrade = "never")
```

Next, install SeuratWrappers in R (as per the [instructions](https://github.com/satijalab/seurat-wrappers)):

```r
# specify the commit number for reproducibility (also cuz later versions require Seurat v5)
remotes::install_github("satijalab/seurat-wrappers", ref = "d28512f", upgrade = "never")
```

#### `CytoTRACE` installation

[CytoTRACE](https://cytotrace.stanford.edu/) was also installed in this environment.

```bash
# download cytotrace
$ wget https://cytotrace.stanford.edu/CytoTRACE_0.3.3.tar.gz
```

```r
# install in R
devtools::install_local("./CytoTRACE_0.3.3.tar.gz", upgrade = "never")

# delete file
file.remove("./CytoTRACE_0.3.3.tar.gz")
```

### `pyscenic`

This environment was created for running pySCENIC (originally created by Shraddha). It contains JupyterLab to run the pySCENIC Jupyter Notebook `.ipynb` files.

### `scrnaseq_env`

This environment was created for running dataset integrations as well as most downstream analysis

Note that for the `scProportionTest` the actual package must be installed separately. First, activate the `scrnaseq_env` environment:

```bash
# activate the scrnaseq_env environment
$ conda activate scrnaseq_env
```

Then install `scProportionTest` in R (as per the [instructions](https://github.com/rpolicastro/scProportionTest/blob/master/README.md) on the scProportionTest Github repo README):

```r
# specify the commit number for reproducibility
remotes::install_github("rpolicastro/scProportionTest", ref = "37a0490")
```

### `clustree`

This environment was created determine the evolution of unsupervised clusters using the `clustree` package.

Also note that it was adapted to conduct pathway analysis as well using the `gprofiler2` package