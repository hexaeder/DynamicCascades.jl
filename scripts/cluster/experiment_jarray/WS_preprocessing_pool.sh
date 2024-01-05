#!/bin/bash

#SBATCH --qos=short
#SBATCH --time=0-04:00:00
#SBATCH --job-name=test_preprocessing_pool
#SBATCH --output=%x-%j-%N.out
#SBATCH --error=%x-%j-%N.err

~/julia-1.8.4/bin/julia WS_preprocessing.jl
