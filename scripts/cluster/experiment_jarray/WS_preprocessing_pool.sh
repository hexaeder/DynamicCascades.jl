#!/bin/bash

#SBATCH --qos=short
#SBATCH --time=1-0:00:00
#SBATCH --job-name=test_params_pool
#SBATCH --output=%x-%j-%N.out
#SBATCH --error=%x-%j-%N.err

~/julia-1.8.4/bin/julia WS_preprocessing.jl
