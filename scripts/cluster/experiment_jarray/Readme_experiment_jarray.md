


** Structure job array framework.
(The framework is the same for WS and RTS, suffixes: HPC for running on PIK-Cluser, pool: HU-Pool)
 - Use Julia 1.8.2 on HPC
   - `cd /home/brandner/tmpjulia/`
   - wget https://julialang-s3.julialang.org/bin/linux/x64/1.8/julia-1.8.2-linux-x86_64.tar.gz
   - tar -xvzf julia-1.8.2-linux-x86_64.tar.gz
   - in submit.sh do `/home/brandner/tmpjulia/julia-1.8.2/bin/julia script.jl   
 - Do `] instantiate` manually on cluster. Doing this through the script, it didn't work.
 - `WS_preprocessing.jl` creates config file and initially stable networks (for WS)
 - `WS_master_experiment.sh` (execute as `./WS_master_experiment.sh`, without `sbatch`)
    creates job array for each inertia value. Slurm parameters need to be adapted
    for each inertia value.
 - `WS_master_experiment_complement_sims.sh`
    (execute as `./WS_master_experiment_complement_sims.sh`, without `sbatch`)
    creates job array for each inertia value and also for every frequency bound.
    This was used when adding missing freuquency bounds to an already existing
    simulation data avoiding to restructure to postprocessing script. In some cases
    its easier to copy around new simulation data by hand. Slurm parameters need
    to be adapted
    for each inertia value and the new indices for the frequency bounds to be simulated
    must be indicated.
 - `WS_job_array_HPC_for_master.sh` invoced by `WS_master_experiment.sh`, runs
   `WS_job.jl`
 - `WS_job.jl` script that is executed for every job.
 - `sacct_postprocessing.sh` creates .txt with slurm informations on the finished jobs (COMPLETED/FAILED, CPUTime, etc.)
 - `WS_postprocessing.jl` Does postprocessing statistics (averages and standard error) and plots
 - -----------------
 - `WS_job_array_HPC.sh` for using job array without complicated `WS_master_experiment.sh`
