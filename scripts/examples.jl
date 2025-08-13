"""
Simple example on how to import WS network, simulate a dynamic cascade and then interactively
plot power flow trajectories of the lines and frequency trajectories of the nodes plus network plot 
with time resolved power flows and frequencies.
"""

include(abspath(@__DIR__, "..", "scripts/helpers_jarray.jl"))


# import Watts-Strogath network
network = import_system(:wattsstrogatz; N=100, k=4, β=0.5, graph_seed=4, μ=0, σ=1, distr_seed=4,
                                K=3, α=0.7, M=0.2u"s^2", γ=0.1u"s", τ=1u"s")

# Run simulation (see various kwargs of `simulate`)
sol = simulate(network;
    initial_fail=[27],
    freq_bound = 0.42);

# ineractive visualization
using NetworkDynamicsInspector
using WGLMakie  

inspect_wrapper(sol; which_trajectories = :failing, tmax=6)




#= Example for rerunning simulations that were done with job array framework:
# NOTE for this to work in `src/DynamicCascades.jl` change `RESULTS_DIR` to the 
directory where the file with the experiment data in it.=#
exp_name_date = "WS_k=4_exp04_vary_I_only_lines_and_nodes_PIK_HPC_K_=3,N_G=32_20250321_171511.976"
initial_fail = 15
task_id = 229

sol = simulate(exp_name_date, task_id, initial_fail);