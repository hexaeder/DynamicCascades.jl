#!/bin/bash

#SBATCH --qos=medium
#SBATCH --time=2-00:00:00
#SBATCH --job-name=230724_01_rtsgmlc_lines_vs_inertia
#SBATCH --output=230724_01_rtsgmlc_lines_vs_inertia.out
#SBATCH --nodes=1
#SBATCH --mail-type=begin        # send email when job begins
#SBATCH --mail-type=end          # send email when job ends
#SBATCH --mail-user=brandner@pik-potsdam.de

module purge
module load julia/1.8.2

julia 230724_01_rtsgmlc_lines_vs_inertia.jl
