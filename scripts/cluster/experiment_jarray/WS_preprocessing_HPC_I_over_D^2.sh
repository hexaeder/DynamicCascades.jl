#!/bin/bash

#SBATCH --qos=priority
#SBATCH --time=1-00:00:00
#SBATCH --job-name=preprocessing_master
#SBATCH --account=icone
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.err


module purge
module load julia/1.8.2

julia WS_preprocessing_I_over_D^2.jl
