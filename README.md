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

Software dependencies are defined in /Manifest.toml and /Project.toml that can be used by Pkg (Julia's builtin package manager) for the installation of the software. The typical install time on a "normal" desktop computer is about 5 to 10 minutes. 
The software has been tested on version `v1.0-submission`.

## Basic functionality
The main entrypoint is `simulate` in src/ND_model.jl.

See `scripts/examples.jl` for how to
 - simulate a cascade in an example network. This takes less than 10 seconds on a "normal" desktop computer. Expected output (up to numerical errors):
```
[ Info: Import system Watts-Strogatz
Find steady state...
Shutdown line 27 at t=0.1
Line 53 tripped at t=4.329652807288284
Line 52 tripped at t=4.547521902036494
Vertex 16 tripped at t=4.598866846813731
Line 60 tripped at t=5.34374202441646
Line 62 tripped at t=5.589868915515144
Vertex 21 tripped at t=5.705406726542723
Terminated on steady state at 65.54382770769196
```
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
 - use the scripts in `scripts/paper_plots/`
 - if you did not download this from Zenodo, download `simulation_results` from [here](https://doi.org/10.5281/zenodo.20341268) and locate them according to the [Structure of the code](#structure-of-the-code)
 - Figure 3 in the paper was created using the
   [legacy branch](https://github.com/hexaeder/DynamicCascades.jl/tree/mwe_old_ND_maybe_plots) that needs an older Julia version and using [Inkscape](https://inkscape.org/). Detailed instructions how to reproduce this figure will be made available upon request.

