#!/bin/bash

#SBATCH --qos=priority
#SBATCH --time=1-00:00:00
#SBATCH --job-name=230728_01_rtsgmlc_lines_vs_inertia_test
#SBATCH --output=230728_01_rtsgmlc_lines_vs_inertia_test.out
#SBATCH --nodes=1
#SBATCH --mail-type=begin        # send email when job begins
#SBATCH --mail-type=end          # send email when job ends
#SBATCH --mail-user=brandner@pik-potsdam.de

module purge
module load julia/1.8.2

julia 230728_01_rtsgmlc_lines_vs_inertia_test.jl
