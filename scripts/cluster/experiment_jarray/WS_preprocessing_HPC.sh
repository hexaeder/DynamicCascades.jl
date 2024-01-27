#!/bin/bash

#SBATCH --qos=priority
#SBATCH --time=0-02:00:00
#SBATCH --job-name=test_preprocessing_master
#SBATCH --account=icone
#SBATCH --output=%x-%j-%N.out
#SBATCH --error=%x-%j-%N.err


module purge
module load julia/1.8.2

julia WS_preprocessing.jl
