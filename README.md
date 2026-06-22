<!-- [![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://wuerfel.io/DynamicCascades.jl/dev/) -->

# DynamicCascades.jl
Exploring dynamic cascades in power systems featuring dynamic line and node failures.

For instructions on how to use the code and reproduce the results of this [paper](https://arxiv.org/abs/2606.20060) see [README.md](https://github.com/hexaeder/DynamicCascades.jl/blob/submission/README.md) in the branch [submission](https://github.com/hexaeder/DynamicCascades.jl/tree/submission).

## Legacy branch
This repository is based on [NetworkDynamics.jl](https://github.com/JuliaDynamics/NetworkDynamics.jl) (ND.jl). NetworkDynamics.jl 
was updated with breaking changes and DynamicCascades.jl was ported to be compatible with the newer ND.jl version. However, only
core functionality was ported. Features that have not been ported are still available on the legacy branch 
[mwe_old_ND_maybe_plots](https://github.com/hexaeder/DynamicCascades.jl/tree/mwe_old_ND_maybe_plots) (branch name: this was initially 
a branch for testing the new ND.jl using a MWE) and are listed below.


### (available) Functionality that has not been ported from mwe_old_ND_maybe_plots
  - asymmetric frequency boundaries
  - different initial perturbation: power perturbation
  - dead band
  - static failures
  - load node failures
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


## Note
 - Don't use steady states from old ND.jl for version with new ND (different order of internal states)
 - Different number of BLAS threads may lead to different results, see scripts/helpers_jarray.jl