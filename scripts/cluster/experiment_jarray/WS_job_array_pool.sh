#!/bin/bash

#SBATCH --qos=short
#SBATCH --time=1-00:00:00
#SBATCH --job-name=test_params_K=9
#SBATCH --output=%x-%A_%a-%N.out
#SBATCH --error=%x-%A_%a-%N.err
#SBATCH --chdir=/users/stud/brandner/MA/repos/MA_data/results_NB/WS_testrun_params_K=6_pool_N_G=2_20240106_021205.818
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --array=1-540

#SBATCH --mail-type=begin        # send email when job begins
#SBATCH --mail-type=end          # send email when job ends
#SBATCH --mail-type=FAIL         # send email when job fails
#SBATCH --mail-user=brandner@pik-potsdam.de

exp_name_date="WS_testrun_params_K=6_pool_N_G=2_20240106_021205.818"

cd /users/stud/brandner/MA/repos/DynamicCascades.jl/scripts/cluster/experiment_jarray

~/julia-1.8.4/bin/julia WS_job.jl $SLURM_ARRAY_TASK_ID $exp_name_date

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
