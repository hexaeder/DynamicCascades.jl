#!/bin/bash

#SBATCH --qos=short
#SBATCH --time=1-00:00:00
#SBATCH --job-name=test_preprocessing
#SBATCH --account=icone
#SBATCH --output=bla.out




echo "Check whether new commits need to be pulled"
echo ""
which git
# /usr/bin/git pull


module purge
module load julia/1.8.2

julia WS_preprocessing.jl
