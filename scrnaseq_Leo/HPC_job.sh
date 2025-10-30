#!/bin/bash
#$ -P pailab
#$ -V
#$ -N liger_sctransform
#$ -cwd
#$ -pe smp 1
#$ -l h_vmem=200G,h_rt=360000
#$ -m beas 
#$ -M llau@oicr.on.ca

# Activate Conda environment
module load conda/24.1.2

# Series env
source /.mounts/labs/pailab/modulator/sw/Ubuntu20.04/conda-24.1.2/bin/activate scrnaseq_env

# Run your script
echo "Starting SingleR"
Rscript utilities/singleR_tumour_predict.R