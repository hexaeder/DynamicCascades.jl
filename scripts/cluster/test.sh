#!/bin/bash

#SBATCH --qos=priority
#SBATCH --time=00:10:00
#SBATCH --job-name=test
#SBATCH --output=test.out
#SBATCH --nodes=1
#SBATCH --mail-type=begin        # send email when job begins
#SBATCH --mail-type=end          # send email when job ends
#SBATCH --mail-user=brandner@pik-potsdam.de

module purge
module load julia/1.8.2

julia test.jl
