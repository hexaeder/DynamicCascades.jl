#!/bin/bash

############################## Parameters to be chosen #########################
name=WS_k=4_exp02_
# inertia_values = [0.2, 0.5, 1.0, 3.0, 5.0, 7.5, 10.0, 20.0, 30.0]
qos_array=(short short short short short short short medium medium)
times_array=(0-05:00:00 0-07:00:00 0-09:00:00 0-11:00:00 0-15:00:00 0-20:00:00 0-22:00:00 1-12:00:00 2-00:00:00)
cpus_array=(1 1 1 1 1 1 1 2 2)

############################## Preprocessing ###################################
echo "------------------------------------------------------------"
echo "Running experiment: $name"
echo "------------------------------------------------------------"

# `cut -f 4 -d' '` extracts the job ID from the output of the `sbatch` command.
# The job ID is stored in the variable PREPROC.
PREPROC=$(sbatch WS_preprocessing_HPC.sh | cut -f 4 -d' ')
echo "SLURM JOB ID Preprocessing: $PREPROC"

sleeptime=400
echo "Sleeping $sleeptime seconds until variables for Slurm are assigned."
sleep $sleeptime

# Read the CSV file into a variable
sbatch_dict=$(cat "/home/brandner/DynamicCascades.jl/scripts/cluster/experiment_jarray/sbatch_dict_$name.csv")
# Extract index value as an array that may be passed to `--array`
exp_name_date=$(echo "$sbatch_dict" | grep 'exp_name_date' | cut -d ',' -f 2- | tr -d '[]' | sed 's/"//g')
job_array_length=$(echo "$sbatch_dict" | grep 'job_array_length' | cut -d ',' -f 2- | tr -d '[]' | sed 's/"//g')
N_inertia=$(echo "$sbatch_dict" | grep 'N_inertia' | cut -d ',' -f 2- | tr -d '[]' | sed 's/"//g')

workdir=/home/brandner/MA_data/results_NB/$exp_name_date/output
sacct_dir=/home/brandner/MA_data/results_NB/$exp_name_date

# NOTE bash starts at index 0.
for job_array_index in $(seq 1 1 $N_inertia); do
    echo "Job array index $job_array_index"
    qos=${qos_array[($job_array_index-1)]}
    time=${times_array[($job_array_index-1)]}
    job_name="$name,idx=$job_array_index"
    cpus=${cpus_array[($job_array_index-1)]}
    ############################## Core Simulation #############################
    MODEL_JOBARRAY=$(sbatch --depend=afterany:$PREPROC --qos=$qos --time=$time --job-name=$job_name --workdir=$workdir --cpus-per-task=$cpus --array=1-$job_array_length WS_job_array_HPC_for_master.sh $exp_name_date $job_array_index | cut -f 4 -d' ')
    echo "SLURM JOB ID JOBARRAY $job_array_index: $MODEL_JOBARRAY"
    ############################## Postprocessing ##############################
    POSTPROC=$(sbatch --depend=afterany:$MODEL_JOBARRAY --job-name="sacct_infos_$job_name" --workdir=$sacct_dir sacct_postprocessing.sh $MODEL_JOBARRAY $sacct_dir $job_array_index | cut -f 4 -d' ')
    echo "SLURM JOB ID Postprocessing $job_array_index: $POSTPROC"
done

echo "Files will be saved in this folder: $exp_name_date"
exit 0

# Explanations

# TODO adapt explanations
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