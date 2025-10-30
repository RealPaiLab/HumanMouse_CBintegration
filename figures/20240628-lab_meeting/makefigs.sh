#!/bin/bash

conda run -n scrnaseq_env Rscript ./makefigs_seurat.R

conda run -n monocle3 Rscript ./makefigs_monocle3.R
