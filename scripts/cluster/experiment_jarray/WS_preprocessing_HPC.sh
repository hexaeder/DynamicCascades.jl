#!/bin/bash

#SBATCH --qos=priority
#SBATCH --time=1-00:00:00
#SBATCH --job-name=WS_k=4_exp05_2_I_over_Dsq_nodes_
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.err


# module purge
# module load julia/1.8.2

# julia WS_preprocessing_I_over_Dsq.jl

/home/brandner/tmpjulia/julia-1.8.4/bin/julia WS_preprocessing_I_over_Dsq_nodes.jl
