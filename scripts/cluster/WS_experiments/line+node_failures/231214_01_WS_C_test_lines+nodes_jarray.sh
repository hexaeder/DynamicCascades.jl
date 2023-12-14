#!/bin/bash

#SBATCH --qos=priority
#SBATCH --time=1-00:00:00
#SBATCH --job-name=231214_01_WS_C_test_lines+nodes_jarray
#SBATCH --output=231214_01_WS_C_test_lines+nodes_jarray.out
#SBATCH --nodes=1
#SBATCH --mail-type=begin        # send email when job begins
#SBATCH --mail-type=end          # send email when job ends
#SBATCH --mail-user=brandner@pik-potsdam.de

module purge
module load julia/1.8.2

julia 231214_01_WS_C_test_lines+nodes_jarray.jl
