#!/bin/bash

#SBATCH --account=icone
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
/home/brandner/tmpjulia/julia-1.8.2/bin/julia WS_job.jl $SLURM_ARRAY_TASK_ID $exp_name_date $job_array_index # $freq_bound_index

# # Specify the path to the config file
# config=231219_WS_test_config.csv

# # extract parameters for the current $SLURM_ARRAY_TASK_ID
# # TODO: avoid code duplification
# beta=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID -F "," '$1==ArrayTaskID {print $2}' $config)
# inertia_values=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID -F "," '$1==ArrayTaskID {print $3}' $config)
# freq_bounds=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID -F "," '$1==ArrayTaskID {print $4}' $config)
# trip_lines=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID -F "," '$1==ArrayTaskID {print $5}' $config)
# trip_nodes=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID -F "," '$1==ArrayTaskID {print $6}' $config)
# k=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID -F "," '$1==ArrayTaskID {print $7}' $config)
# N_nodes=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID -F "," '$1==ArrayTaskID {print $8}' $config)
# K=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID -F "," '$1==ArrayTaskID {print $9}' $config)
# gamma=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID -F "," '$1==ArrayTaskID {print $10}' $config)
# tau=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID -F "," '$1==ArrayTaskID {print $11}' $config)
# alpha=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID -F "," '$1==ArrayTaskID {print $12}' $config)
# init_pert=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID -F "," '$1==ArrayTaskID {print $13}' $config)
# sigma=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID -F "," '$1==ArrayTaskID {print $14}' $config)
# mu=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID -F "," '$1==ArrayTaskID {print $15}' $config)

# For testing: print output.txt the current $SLURM_ARRAY_TASK_ID and corresponding parameters
# echo "This is array task ${SLURM_ARRAY_TASK_ID}, with parameters \n
#     beta=${beta},\n
#     inertia=${inertia_values},\n
#     frequency bound=${freq_bounds},\n
#     trip_lines=${trip_lines},\n
#     trip_nodes=${trip_nodes},\n
#     k=${k},\n
#     N_nodes=${N_nodes},\n
#     K=${K},\n
#     gamma=${gamma},\n
#     alpha=${alpha},\n
#     init_pert=${init_pert},\n
#     sigma=${sigma},\n
#     mu=${mu}." >> output.txt
