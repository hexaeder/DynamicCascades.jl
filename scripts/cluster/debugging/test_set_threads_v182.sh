#!/bin/bash

#SBATCH --qos=priority
#SBATCH --time=24:00:00
#SBATCH --job-name=test_set_threads_v182
#SBATCH --output=test_set_threads_v182.out
#SBATCH --nodes=1

#SBATCH --mail-type=begin        # send email when job begins
#SBATCH --mail-type=end          # send email when job ends
#SBATCH --mail-user=brandner@pik-potsdam.de

module purge
module load julia/1.8.2

julia test_set_threads_v182.jl
