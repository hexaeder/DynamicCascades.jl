#!/bin/bash

#SBATCH --qos=priority
#SBATCH --time=0-01:00:00
#SBATCH --account=icone
#SBATCH --output=%x-%A_%a-%N.out
#SBATCH --error=%x-%A_%a-%N.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1

SLURM_JOBID=$1
cd $2
job_array_index=$3

sacct --format=JobID,JobName,Submit,Start,End,State,Partition,ReqTRES%30,CPUTime,Elapsed,MaxRSS,NodeList%30 --units=M -j $SLURM_JOBID > "sacct_infos_idx=${job_array_index},${SLURM_JOBID}.txt"
