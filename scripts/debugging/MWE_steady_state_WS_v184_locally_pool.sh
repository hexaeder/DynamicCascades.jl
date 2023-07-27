#!/bin/bash

#SBATCH --qos=priority
#SBATCH --time=24:00:00
#SBATCH --job-name=MWE_steady_state_WS_v184_locally_pool
#SBATCH --output=MWE_steady_state_WS_v184_locally_pool.out
#SBATCH --nodes=1

#SBATCH --mail-type=begin        # send email when job begins
#SBATCH --mail-type=end          # send email when job ends
#SBATCH --mail-user=brandner@pik-potsdam.de

srun ~/tmpjulia/julia-1.8.4/bin/julia -t 1 MWE_steady_state_WS_v184_locally_pool.jl
