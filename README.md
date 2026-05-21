# DynamicCascades.jl
Exploring dynamic cascades in power systems featuring dynamic line and node failures.


## Structure of the code
submission/
 - data/ network data and parameters for the RTS-GMLC test grid (RTS)
 - generated_figs/
 - scripts/
    - paper_plots/ scripts to generate the figures in the paper, they are saved in the directory `generated figs`
    - RTS/ scripts to generate the results for RTS contained in the directory `simulation_results`
    - WS/ scripts to generate the results for Watts-Strogatz networks (WS) contained in the directory `simulation_results`
    - examples.jl simple examples how to run the code
 - simulation_results/ 
 - src/ core source code defining the dynamic cascading failure model


## Basic functionality
The main entrypoint is `simulate` in src/ND_model.jl.

See `scripts/examples.jl` for how to
 - simulate a cascade in an example network
 - resimulate a cascade that has been simulated in an experiment using the job array framework (see [below](#reproduction-of-the-results-in-simulation_results)) 
 - use the interactive plotting utility [NetworkDynamicsInspector.jl](https://github.com/JuliaDynamics/NetworkDynamics.jl/tree/main/NetworkDynamicsInspector) 


## Reproduction of the results in simulation_results
This uses a job array framework for High-performance computing (HPC). The framework is the same for WS and RTS.

 - `WS_preprocessing.jl` creates config file and initially stable networks (for WS).
    Adapt the parameters in this script for each new experiment.
 - `WS_master_experiment.sh` (execute as `./WS_master_experiment.sh`, without `sbatch`)
    creates separate job array for each inertia value to reduce cluster queue time. Slurm parameters need
    to be adapted for each inertia value.
 - `WS_job_array_HPC_for_master.sh` invoked by `WS_master_experiment.sh`, runs `WS_job.jl`
 - `WS_job.jl` script that is executed for every job.
 - `sacct_postprocessing.sh` creates .txt with slurm information on the finished jobs (COMPLETED/FAILED, CPUTime, etc.)


## Reproduction of the figures in the paper
 - create the directory `generated_figs` (see above)
<!-- - if you did not download this from Zenodo, download `simulation_results` from [here](TODO) and locate them according to the [Structure of the code](#structure-of-the-code) -->
 - Figure 3 in the paper was created using the
   [legacy branch](https://github.com/hexaeder/DynamicCascades.jl/tree/mwe_old_ND_maybe_plots) that needs an older Julia version and using [Inkscape](https://inkscape.org/). Detailed instructions how to reproduce this figure will be made available upon request.

