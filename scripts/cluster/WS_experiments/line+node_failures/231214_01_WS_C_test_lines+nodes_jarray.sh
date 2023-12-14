#!/bin/bash
#SBATCH --job-name=myarrayjob
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --array=1-2


# Specify the path to the config file
config=config.txt

# extract parameters for the current $SLURM_ARRAY_TASK_ID
sample=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $2}' $config)


sex=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $3}' $config)

# Print to a file a message that includes the current $SLURM_ARRAY_TASK_ID, the same name, and the sex of the sample
echo "This is array task ${SLURM_ARRAY_TASK_ID}, the sample name is ${sample} and the sex is ${sex}." >> output.txt

module purge
module load julia/1.8.2

julia test_jarray.jl $sex
# julia 231214_01_WS_C_test_lines+nodes_jarray.jl
