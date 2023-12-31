#!/bin/bash
#SBATCH --qos=short
#SBATCH --time=1-00:00:00
#SBATCH --job-name=test_preprocessing
#SBATCH --account=icone
#SBATCH --output=test-%x-%j-%N.out
#SBATCH --error=test-%x-%j-%N.err

module purge
module load julia/1.8.2

julia WS_preprocessing.jl
