include(abspath(@__DIR__, "..", "helpers_jarray.jl"))
include(abspath(@__DIR__, "WS_trajectories_new_ND_single_model_port.jl"))

using DynamicCascades
using NetworkDynamics
using Serialization

# Load the data
exp_name_date = "WS_k=4_exp04_vary_I_only_lines_and_nodes_PIK_HPC_K_=3,N_G=32_20250321_171511.976"
# exp_name_date = "WS_k=4_exp07_3_vary_D_lines_and_nodes_PIK_HPC_K_=3,N_G=32_20250326_231211.092"
# exp_name_date = "WS_k=4_exp03_2_vary_I_only_nodes_PIK_HPC_K_=3,N_G=32_20250326_231746.273"
exp_data_dir = joinpath(RESULTS_DIR, exp_name_date)
df_config = DataFrame(CSV.File(joinpath(exp_data_dir, "config.csv")));

###
### run multiple simulations
###
# task_id_array = [2307, 2309, 2312, 2313] # inertia variation
# task_id_array = [908, 912, 916, 920, 924, 928, 932, 936] # damping variation
task_id_array = [479, 481, 483, 485]
initial_fail_array = [95]
for task_id in task_id_array
    for initial_fail in initial_fail_array
        println("task_id=$task_id, initial_fail=$initial_fail")
        sol, nw = simulate_new_ND(exp_data_dir, task_id, initial_fail;
            solverargs = (;dtmax=0.01),
            verbose = true);
        Serialization.serialize(joinpath(exp_data_dir, "general_investigations", "task_id=$task_id,initial_fail=$initial_fail.sol"), sol)
    end
end

###
### run single solution
###
task_id = 2310
initial_fail = 95
sol = simulate_new_ND(exp_data_dir, task_id, initial_fail;
    solverargs = (;dtmax=0.01),
    verbose = true);

Serialization.serialize(joinpath(exp_data_dir, "general_investigations", "task_id=$task_id,initial_fail=$initial_fail.sol"), sol)

###
### plot 
###
using Colors
using CairoMakie
plot_dir = simulation_dir = joinpath(exp_data_dir, "general_investigations")

# single solution
plot_simulation(simulation_dir, plot_dir, df_config, task_id, initial_fail)

# multiple solutions
# task_id_array = [2307, 2309, 2312, 2313] # inertia variation
# task_id_array = [908, 912, 916, 920, 924, 928, 932, 936] # damping variation
task_id_array = [479, 481, 483, 485] # nodes only inertia variation
for task_id in task_id_array
    println("task_id=$task_id, initial_fail=$initial_fail")
    plot_simulation(simulation_dir, plot_dir, df_config, task_id, initial_fail)
end


###
### inspect solutions
###
using NetworkDynamicsInspector
using WGLMakie # for inspector
# using Bonito # for using plot pane and memorizing plots
# Bonito.set_cleanup_time!(720)
task_id = 1446
initial_fail = 15

sol = Serialization.deserialize(joinpath(exp_data_dir, "general_investigations", "task_id=$task_id,initial_fail=$initial_fail.sol"));
nw = Serialization.deserialize(joinpath(exp_data_dir, "general_investigations", "task_id=$task_id,initial_fail=$initial_fail.nw"));
inspect(sol; restart=true, reset=true)
set_sol!(sol) # optional if after inspect(sol)
set_graphplot!(; nstate=[:ω], estate=[:S], nstate_rel=false, estate_rel=false)
vindices = [i for i in 1:nv(nw) if sol(sol.t[end], idxs=vidxs(i, :node_swing_stat))[1] == 0]
eindices = [i for i in 1:ne(nw) if sol(sol.t[end], idxs=eidxs(i, :line_stat))[1] == 0]
define_timeseries!([
    (; selcomp=[VIndex(i) for i in vindices], states=[:ω, :ωmax], rel=false),
    (; selcomp=[EIndex(i) for i in eindices], states=[:S, :rating], rel=false),
])

ensemble_element=9,I=7.5,D=1,f_b=0.045,task_id=2310,initial_fail=95