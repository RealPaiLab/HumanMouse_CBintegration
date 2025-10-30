#!/bin/bash
#$ -P pailab
#$ -V
#$ -N human_ubc_pyscenic
#$ -cwd
#$ -pe smp 24
#$ -l h_vmem=8G,h_rt=1:0:0:0
#$ -m beas
#$ -M icheong@oicr.on.ca

outRoot=/.mounts/labs/pailab/private/projects/MB_multiome/output/PySCENIC

# print date yymdd
dt=$(date +%y%m%d)
outDir=${outRoot}/${dt}
mkdir -p $outDir

# activate conda environment
module load conda/24.1.2
source activate base
conda activate scrnaseq_env

# extract raw counts, metadata, and UMAP embeddings from Seurat object

Rscript seurat_to_scenic.R 


#### run pySCENIC
###
###conda activate pyscenic
###
###python /.mounts/labs/pailab/private/icheong/CBL_scRNAseq/HumanMouseUBC/integrated_human_ubc/03_pyscenic.py \
###    --counts_csv $out_dir/human_ubcs.rna_counts.csv \
###    --metadata_csv $out_dir/human_ubcs.metadata.csv \
###    --umap_csv $out_dir/human_ubcs.umap_embeddings.csv \
###    --cluster_column SCT_snn_res.0.5 \
###    --num_threads 24 \
###    --dataset_id human_ubcs \
###    --out_dir $out_dir \
###    &> run_pyscenic.log
###
###conda deactivate
