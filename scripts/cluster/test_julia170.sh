#!/bin/bash

#SBATCH --qos=priority
#SBATCH --time=24:00:00
#SBATCH --job-name=test_julia170
#SBATCH --output=test_julia170.out
#SBATCH --nodes=1
#SBATCH --mail-type=begin        # send email when job begins
#SBATCH --mail-type=end          # send email when job ends
#SBATCH --mail-user=brandner@pik-potsdam.de

module purge
module load julia/1.7.0

julia test_julia170.jl
