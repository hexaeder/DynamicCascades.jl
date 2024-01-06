#!/bin/bash

#SBATCH --qos=short
#SBATCH --time=0-02:00:00
#SBATCH --job-name=test_params_pool_K=9
#SBATCH --output=%x-%j-%N.out
#SBATCH --error=%x-%j-%N.err

~/julia-1.8.4/bin/julia WS_preprocessing.jl
