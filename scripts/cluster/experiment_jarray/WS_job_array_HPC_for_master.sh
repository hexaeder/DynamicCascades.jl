#!/bin/bash

#SBATCH --output=%x-%A_%a-%N.out
#SBATCH --error=%x-%A_%a-%N.err
#SBATCH --ntasks=1

#SBATCH --mail-type=begin        # send email when job begins
#SBATCH --mail-type=end          # send email when job ends
#SBATCH --mail-type=FAIL         # send email when job fails
#SBATCH --mail-user=brandner@pik-potsdam.de

exp_name_date=$1 # getting variable from master script
job_array_index=$2

# This is only for `WS_master_experiment_complement_sims.sh`
# freq_bound_index=$3

cd /home/brandner/DynamicCascades.jl/scripts/cluster/experiment_jarray

# module purge
# module load julia/1.8.2

# `$freq_bound_index` is only for `WS_master_experiment_complement_sims.sh`
# julia WS_job.jl $SLURM_ARRAY_TASK_ID $exp_name_date $job_array_index # $freq_bound_index
/home/brandner/tmpjulia/julia-1.8.4/bin/julia WS_job.jl $SLURM_ARRAY_TASK_ID $exp_name_date $job_array_index # $freq_bound_index
