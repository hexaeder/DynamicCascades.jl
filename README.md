<!-- [![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://wuerfel.io/DynamicCascades.jl/dev/) -->

# DynamicCascades.jl
Exploring dynamic cascades in power systems featuring dynamic line and node failures.

## Basic functionality
The main entrypoint is `simulate` in src/ND_model.jl.

See `scripts/examples.jl` for how to
 - simulate a cascade in an example network
 - resimulate a cascade that has been simulated in an experiment using the job array framework (see below) 
 - use the interactive plotting utility [NetworkDynamicsInspector](https://github.com/JuliaDynamics/NetworkDynamics.jl/tree/main/NetworkDynamicsInspector) 

##  Job array framework for HPCs
The framework is the same for WS and RTS.

<!--
 - For reproduction use Julia 1.8.4 (for WS) and Julia v1.11.0 (for RTS) on HPC 
   - Example for Julia 1.8.2 
     - `cd /home/brandner/tmpjulia/`
     - wget https://julialang-s3.julialang.org/bin/linux/x64/1.8/julia-1.8.2-linux-x86_64.tar.gz
     - tar -xvzf julia-1.8.2-linux-x86_64.tar.gz
     - in submit.sh do `/home/brandner/tmpjulia/julia-1.8.2/bin/julia script.jl   
     - Do `] instantiate` manually on cluster. Doing this through the script, it didn't work. 
-->

 - `WS_preprocessing.jl` creates config file and initially stable networks (for WS).
    Adapt the parameters in this script for each new experiment. No separate script for each 
    experiment needed as the script that was used for the simulations is saved with the 
    simulation data when executing `WS_master_experiment.sh` and `RTS_master_experiment.sh`.
 - `WS_master_experiment.sh` (execute as `./WS_master_experiment.sh`, without `sbatch`)
    creates separate job array for each inertia value to reduce cluster queue time. Slurm parameters need
    to be adapted for each inertia value.
 - `WS_job_array_HPC_for_master.sh` invoced by `WS_master_experiment.sh`, runs `WS_job.jl`
 - `WS_job.jl` script that is executed for every job.
 - `sacct_postprocessing.sh` creates .txt with slurm informations on the finished jobs (COMPLETED/FAILED, CPUTime, etc.)

## Legacy branch
This repository is based on [NetworkDynamics.jl](https://github.com/JuliaDynamics/NetworkDynamics.jl) (ND.jl). NetworkDynamics.jl 
was updated with breaking changes and DynamicCascades.jl was ported to be compatible with the newer ND.jl version. However, only
core functionality was ported. Features that have not been ported are still available on the legacy branch 
[mwe_old_ND_maybe_plots](https://github.com/hexaeder/DynamicCascades.jl/tree/mwe_old_ND_maybe_plots) (branch name: this was initially 
a branch for testing the new ND.jl using a MWE) and are listed below.


### (available) Functionality that has not been portet from mwe_old_ND_maybe_plots
  - assymmetric frequency boundaries
  - different initial perturbations: power perturbation, node failure (?)
  - dead band
  - static failures
  - load node falures
  - check active vs apparent power flow in line CB
  - Further networks
    - import_isolator_toymodel.jl
    - import_kaiser2020.jl
    - import_rts96.jl
    - import_schaefer2018.jl
    - import_slack_gen.jl
    - import_square.jl
    - import_toymodel.jl
  - Plotting functionalities in inspect_solution.jl
    - tbd: mention snapshot plot 


<!-- Some examples using deprecated code can be found in [here](https://wuerfel.io/DynamicCascades.jl/dev/). -->