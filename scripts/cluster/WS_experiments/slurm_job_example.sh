#!/bin/bash
#SBATCH --qos=medium
#SBATCH --time=3-00:00:00
#SBATCH --job-name=test_job
#SBATCH --account=icone
#SBATCH --output=test-%x-%j-%N.out
#SBATCH --error=test-%x-%j-%N.err
#SBATCH --workdir=/p/tmp/<your username>/<some directory>

#SBATCH --mail-type=begin        # send email when job begins
#SBATCH --mail-type=end          # send email when job ends
#SBATCH --mail-user=brandner@pik-potsdam.de

echo "------------------------------------------------------------"
echo "SLURM JOB ID: $SLURM_JOBID"
echo "$SLURM_NTASKS tasks"
echo "------------------------------------------------------------"

# TODO maybe use workdir as used in paper simulations

hostname
