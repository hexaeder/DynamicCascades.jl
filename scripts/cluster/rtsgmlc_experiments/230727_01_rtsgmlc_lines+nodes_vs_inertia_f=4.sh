#!/bin/bash

#SBATCH --qos=medium
#SBATCH --time=3-00:00:00
#SBATCH --job-name=230727_01_rtsgmlc_lines+nodes_vs_inertia_f=4
#SBATCH --output=230727_01_rtsgmlc_lines+nodes_vs_inertia_f=4.out
#SBATCH --nodes=1
#SBATCH --mail-type=begin        # send email when job begins
#SBATCH --mail-type=end          # send email when job ends
#SBATCH --mail-user=brandner@pik-potsdam.de

module purge
module load julia/1.8.2

julia 230727_01_rtsgmlc_lines+nodes_vs_inertia_f=4.jl
