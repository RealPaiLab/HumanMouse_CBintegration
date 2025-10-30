#!/bin/bash

#$ -P pailab
#$ -V
#$ -N integ_dataset_CCA
#$ -cwd
#$ -pe smp 7
#$ -l h_vmem=50G,h_rt=300000
#$ -m beas 
#$ -M llau@oicr.on.ca

# Activate Conda environment
module load conda/24.1.2

# Parallel env
source /.mounts/labs/pailab/modulator/sw/Ubuntu20.04/conda-24.1.2/bin/activate scrnaseq_parallel_env

# Setting destination location
destination_date=$(date +'%Y%m%d')
destination_dir="/.mounts/labs/pailab/private/llau/results/integrated/${destination_date}"

# Run your script
echo "Starting integration"
mkdir "$destination_dir"
cp ./HPC_job_parallel.sh "${destination_dir}"
Rscript HPC_integrate.R -datasets Aldinger_RL_human Sepp_RL_human Luo_RL_human Vladoiu_RL_mouse Sepp_RL_mouse --CCA_int --run_series 
mv ./integ_dataset_CCA* "${destination_dir}"
