#!/bin/bash

#SBATCH --qos=short
#SBATCH --time=1-00:00:00
#SBATCH --job-name=test_preprocessing
#SBATCH --account=icone
#SBATCH --output=%x-%j-%N.out


# For having seperate error file #SBATCH --error=%x-%j-%N.err

echo "Check whether new commits need to be pulled"
echo ""
/usr/bin/git/git pull


module purge
module load julia/1.8.2

julia WS_preprocessing.jl
