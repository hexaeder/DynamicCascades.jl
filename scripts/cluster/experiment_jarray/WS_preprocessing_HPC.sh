#!/bin/bash

#SBATCH --qos=priority
#SBATCH --time=0-04:00:00
#SBATCH --job-name=test_preprocessing
#SBATCH --account=icone
#SBATCH --output=%x-%j-%N.out
#SBATCH --error=%x-%j-%N.err


module purge
module load julia/1.8.2

julia WS_preprocessing.jl