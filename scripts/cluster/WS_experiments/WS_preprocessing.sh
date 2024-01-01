#!/bin/bash

#SBATCH --qos=short
#SBATCH --time=1-00:00:00
#SBATCH --job-name=test_preprocessing
#SBATCH --account=icone
#SBATCH --output=%x-%j-%N.out
#SBATCH --error=%x-%j-%N.err


# # TODO This is more complicated as this requires internet.
# echo "Check whether new commits need to be pulled"
# module load git
# git pull



module purge
module load julia/1.8.2

julia WS_preprocessing.jl
