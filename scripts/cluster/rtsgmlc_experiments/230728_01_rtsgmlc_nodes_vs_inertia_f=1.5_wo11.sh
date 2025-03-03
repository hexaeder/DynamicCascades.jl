#!/bin/bash

#SBATCH --qos=medium
#SBATCH --time=5-00:00:00
#SBATCH --job-name=230728_01_rtsgmlc_nodes_vs_inertia_f=1.5_wo11
#SBATCH --output=230728_01_rtsgmlc_nodes_vs_inertia_f=1.5_wo11.out
#SBATCH --nodes=1
#SBATCH --mail-type=begin        # send email when job begins
#SBATCH --mail-type=end          # send email when job ends
#SBATCH --mail-user=brandner@pik-potsdam.de

module purge
module load julia/1.8.2

julia 230728_01_rtsgmlc_nodes_vs_inertia_f=1.5_wo11.jl
