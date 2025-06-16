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

cd /home/brandner/DynamicCascades.jl/scripts/RTS

julia RTS_job.jl $SLURM_ARRAY_TASK_ID $exp_name_date $job_array_index
