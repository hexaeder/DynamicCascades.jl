#!/bin/bash

#SBATCH --qos=short
#SBATCH --time=24:00:00
#SBATCH --job-name=230722_01_rtsgmlc_failures_vs_inertia
#SBATCH --output=230722_01_rtsgmlc_failures_vs_inertia.out
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16


julia -p 16 230722_01_rtsgmlc_failures_vs_inertia.jl
