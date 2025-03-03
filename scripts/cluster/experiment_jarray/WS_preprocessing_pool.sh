#!/bin/bash

#SBATCH --qos=short
#SBATCH --time=1-00:00:00
#SBATCH --job-name=test_params_pool_K=3
#SBATCH --output=%x-%j-%N.out
#SBATCH --error=%x-%j-%N.err

~/julia-1.8.4/bin/julia WS_preprocessing.jl
