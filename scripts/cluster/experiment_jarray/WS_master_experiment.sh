#!/bin/bash

############################## Preprocessing ###################################
name=WS_k=4_exp01_ # NOTE: use `name` from WS_preprocessing.jl!
echo "------------------------------------------------------------"
echo "Running experiment: $name"
echo "------------------------------------------------------------"

# `cut -f 4 -d' '` extracts the job ID from the output of the `sbatch` command.
# The job ID is stored in the variable PREPROC.
PREPROC=$(sbatch WS_preprocessing_HPC.sh | cut -f 4 -d' ')
echo "SLURM JOB ID Preprocessing: $PREPROC"

echo "Sleeping 200 seconds until variables for Slurm are assigned."
sleep 200

output_file=$(cat "/home/brandner/DynamicCascades.jl/scripts/cluster/experiment_jarray/$name.csv")
# Extract index value as an array that may be passed to `--array`
exp_name_date=$(echo "$output_file" | grep 'exp_name_date' | cut -d ',' -f 2- | tr -d '[]' | sed 's/"//g')


## Read in indices in bash-script
# Read the CSV file into a variable
sbatch_dict=$(cat "/home/brandner/MA_data/results_NB/$exp_name_date/sbatch_dict.csv")

# Extract index value as an array that may be passed to `--array` (explanations see below)
indices_short=$(echo "$sbatch_dict" | grep 'indices_short' | cut -d ',' -f 2- | tr -d '[]' | sed 's/"//g')
indices_long=$(echo "$sbatch_dict" | grep 'indices_long' | cut -d ',' -f 2- | tr -d '[]' | sed 's/"//g')

workdir=/home/brandner/MA_data/results_NB/$exp_name_date/output
sacct_dir=/home/brandner/MA_data/results_NB/$exp_name_date

############################## Core Simulation #################################
# Submit job arrays
# https://stackoverflow.com/questions/77479142/slurm-array-add-variable-in-sbatch-options?noredirect=1&lq=1
MODEL_JOBARRAY1=$(sbatch --depend=afterany:$PREPROC --qos=short --time=1-00:00:00 --job-name=WS_k=4_exp01_short --workdir=$workdir --array="${indices_short}" --cpus-per-task=1 WS_job_array_HPC_for_master.sh $exp_name_date | cut -f 4 -d' ')
echo "SLURM JOB ID JOBARRAY 1: $MODEL_JOBARRAY1"

MODEL_JOBARRAY2=$(sbatch --depend=afterany:$PREPROC --qos=medium --time=2-00:00:00 --job-name=WS_k=4_exp01_medium --workdir=$workdir --array="${indices_long}" --cpus-per-task=2 WS_job_array_HPC_for_master.sh $exp_name_date | cut -f 4 -d' ')
echo "SLURM JOB ID JOBARRAY 2: $MODEL_JOBARRAY2"

echo "Files will be saved in this folder: $exp_name_date"

############################## Postprocessing ##################################
POSTPROC1=$(sbatch --depend=afterany:$MODEL_JOBARRAY1 --workdir=$sacct_dir WS_sacct_postprocessing.sh $MODEL_JOBARRAY1 $sacct_dir | cut -f 4 -d' ')
echo "SLURM JOB ID Postprocessing 1: $POSTPROC1"
POSTPROC2=$(sbatch --depend=afterany:$MODEL_JOBARRAY2 --workdir=$sacct_dir WS_sacct_postprocessing.sh $MODEL_JOBARRAY2 $sacct_dir | cut -f 4 -d' ')
echo "SLURM JOB ID Postprocessing 1: $POSTPROC2"

exit 0



# Explanations

# Explanations for `indices_short=$(echo "$sbatch_dict" | grep 'indices_short' | cut -d ',' -f 2- | tr -d '[]' | sed 's/"//g')`
# `echo "$sbatch_dict"` Echo the content of the variable sbatch_dict to standard output.
# `grep 'indices_short'` Use grep to filter lines containing the string 'indices_short'.
# `cut -d ',' -f 2-` Use cut to split the line into fields using a comma as the delimiter, keeping fields starting from the second field.
# This removes the identifier (like ':indices_short') and retains the values.
# `tr -d '[]'` Use tr to delete characters specified by the set '[]'.
#  This removes square brackets, which might be present around the values in the CSV data.
# `sed 's/"//g'` Use sed for text stream editing. The expression 's/"//g' specifies a substitution command to replace all double quotes with an empty string.
# This is used to remove any double quotes around the values.
# `indices_short=$(...)`: This is command substitution. It assigns the result of the entire command sequence within the parentheses to the variable indices_short.
