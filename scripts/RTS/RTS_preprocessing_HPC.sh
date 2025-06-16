#!/bin/bash

#SBATCH --qos=priority
#SBATCH --time=1-00:00:00
#SBATCH --job-name=RTS_exp04_variation_frequency+inertia_
#SBATCH --account=icone
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.err

julia RTS_preprocessing.jl
