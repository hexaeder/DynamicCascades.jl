#!/bin/bash

# Create array of indices written to a file array_indices.txt in preprocessing script
# Maybe also
    # write times to file and hand them to WS_job_array_HPC.sh
    # write experiment file name

# Read in indices in bash-script

# Read the CSV file into a variable
sbatch_dict=$(cat /home/brandner/MA_data/results_NB/WS_testrun_masterPIK_HPC_K_=3,N_G=2_20240127_104017.493/sbatch_dict.csv)

# Extract the value of :indices_short as an array without using jq
indices_short=$(echo "$sbatch_dict" | grep 'indices_short' | cut -d ',' -f 2- | tr -d '[]')
indices_long=$(echo "$sbatch_dict" | grep 'indices_long' | cut -d ',' -f 2- | tr -d '[]')
echo $indices_short
echo $indices_long

# indices_long=

# https://stackoverflow.com/questions/77479142/slurm-array-add-variable-in-sbatch-options?noredirect=1&lq=1
# Submit job arrays
sbatch --array=${indices_short} WS_job_array_HPC.sh
# sbatch WS_job_array_HPC.sh
