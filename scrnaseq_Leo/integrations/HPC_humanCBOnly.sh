#!/bin/bash

#$ -P pailab
#$ -V
#$ -N integ_humanCB_CCA
#$ -cwd
#$ -pe smp 7
#$ -l h_vmem=50G,h_rt=300000
#$ -m beas 
#$ -M spai@oicr.on.ca


# FromLeo: For the RL integration, you would probably only need like 20 GB max per core so 140 GB total memory. For FC integration i would keep it at 50 GB per core just in case.

# Time considerations: CBL with all the datasets took like 2 days, but with only 2 I would guess would be a couple of hours (6 hours)
# RL should be faster

# Leave h_rt as is.

# SP rerunning integration of just human datasets
# Activate Conda environment
module load conda/24.1.2

# Parallel env
source /.mounts/labs/pailab/modulator/sw/Ubuntu20.04/conda-24.1.2/bin/activate scrnaseq_parallel_env

# Setting destination location
destination_date=$(date +'%Y%m%d')
destination_dir="/.mounts/labs/pailab/private/projects/HumanMouseUBC/results/integrated_HsRLonly/${destination_date}"

# Run your script
echo "Starting integration"
mkdir "$destination_dir"
cp ./HPC_job_parallel.sh "${destination_dir}"
Rscript HPC_integrate.R -datasets Aldinger_full_cerebellum_human Sepp_full_cerebellum_human --CCA_int --run_series 
#mv ./integ_dataset_CCA* "${destination_dir}"
 