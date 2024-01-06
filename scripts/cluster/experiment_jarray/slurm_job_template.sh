#!/bin/bash

#SBATCH --qos=short
#SBATCH --time=1-00:00:00
#SBATCH --job-name=test_job
#SBATCH --account=icone
#SBATCH --output=%x-%j-%N.out
#SBATCH --error=%x-%j-%N.err
#SBATCH --workdir=/p/tmp/<your username>/<some directory>

#SBATCH --mail-type=begin        # send email when job begins
#SBATCH --mail-type=end          # send email when job ends
#SBATCH --mail-type=FAIL         # send email when job fails
#SBATCH --mail-user=brandner@pik-potsdam.de


# # TODO This is more complicated as this requires internet.
# echo "Check whether new commits need to be pulled"
# module load git
# git pull


echo "------------------------------------------------------------"
echo "SLURM JOB ID: $SLURM_JOBID"
echo "$SLURM_NTASKS tasks"
echo "------------------------------------------------------------"

# TODO maybe use workdir as used in paper simulations

hostname

module purge
module load julia/1.8.2

julia test.jl
