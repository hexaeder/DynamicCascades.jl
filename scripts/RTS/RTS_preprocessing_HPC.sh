#!/bin/bash

#SBATCH --qos=priority
#SBATCH --time=1-00:00:00
#SBATCH --job-name=RTS_test_new_ND_
#SBATCH --account=icone
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.err


# module purge
# module load julia/1.8.2

julia RTS_preprocessing.jl
