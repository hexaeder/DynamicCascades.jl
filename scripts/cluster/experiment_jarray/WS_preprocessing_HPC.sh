#!/bin/bash

#SBATCH --qos=priority
#SBATCH --time=1-00:00:00
#SBATCH --job-name=preprocessing_WS_k=4_exp06_1_I_over_Dsq_lines__
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.err


# module purge
# module load julia/1.8.2

# julia WS_preprocessing_I_over_Dsq.jl

/home/brandner/tmpjulia/julia-1.8.4/bin/julia WS_preprocessing_I_over_D_lines.jl
