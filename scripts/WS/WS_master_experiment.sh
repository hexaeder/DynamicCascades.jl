#!/bin/bash

############################## Parameters to be chosen #########################
name=WS_k=4_exp11_vary_I_only_lines_and_nodes_change_to_BH_complement_f_b
# inertia_values = [0.2, 0.5, 1.0, 3.0, 5.0, 7.5, 10.0, 20.0, 30.0]
qos_array=(short short short short short short short short short)
times_array=(0-05:00:00 0-07:00:00 0-09:00:00 0-11:00:00 0-15:00:00 0-20:00:00 0-22:00:00 0-24:00:00 0-24:00:00)
cpus_array=(1 1 1 1 1 1 1 1 1)

############################## Preprocessing ###################################
echo "------------------------------------------------------------"
echo "Running experiment: $name"
echo "------------------------------------------------------------"

# `cut -f 4 -d' '` extracts the job ID from the output of the `sbatch` command.
# The job ID is stored in the variable PREPROC.
PREPROC=$(sbatch WS_preprocessing_HPC.sh | cut -f 4 -d' ')
echo "SLURM JOB ID Preprocessing: $PREPROC"

sleeptime=300
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
cp -r /home/brandner/DynamicCascades.jl $workdir

# NOTE bash starts at index 0.
for job_array_index in $(seq 1 1 $N_inertia); do
    echo "Job array index $job_array_index"
    qos=${qos_array[($job_array_index-1)]}
    time=${times_array[($job_array_index-1)]}
    job_name="$name,idx=$job_array_index"
    cpus=${cpus_array[($job_array_index-1)]}
    ############################## Core Simulation #############################
    MODEL_JOBARRAY=$(sbatch --depend=afterany:$PREPROC --qos=$qos --time=$time --job-name=$job_name --chdir=$workdir --cpus-per-task=$cpus --array=1-$job_array_length WS_job_array_HPC_for_master.sh $exp_name_date $job_array_index | cut -f 4 -d' ')
    echo "SLURM JOB ID JOBARRAY $job_array_index: $MODEL_JOBARRAY"
    ############################## Postprocessing ##############################
    POSTPROC=$(sbatch --depend=afterany:$MODEL_JOBARRAY --job-name="sacct_infos_$job_name" --chdir=$sacct_dir sacct_postprocessing.sh $MODEL_JOBARRAY $sacct_dir $job_array_index | cut -f 4 -d' ') 
    echo "SLURM JOB ID Postprocessing $job_array_index: $POSTPROC"
done

echo "Files will be saved in this folder: $exp_name_date"
exit 0