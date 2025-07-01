#!/bin/bash

#SBATCH --qos=priority
#SBATCH --time=1-00:00:00
#SBATCH --job-name=WS_k=4_exp11_vary_I_only_lines_and_nodes_change_to_BH_complement_f_b
#SBATCH --account=icone
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.err

julia WS_preprocessing.jl
