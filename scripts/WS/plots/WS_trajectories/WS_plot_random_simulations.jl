include(abspath(@__DIR__, "..", "helpers_jarray.jl"))


using DynamicCascades
using NetworkDynamics
using Serialization

# Load the data
exp_name_date = "WS_k=4_exp04_vary_I_only_lines_and_nodes_PIK_HPC_K_=3,N_G=32_20250321_171511.976"
exp_data_dir = joinpath(RESULTS_DIR, exp_name_date)
df_config = DataFrame(CSV.File(joinpath(exp_data_dir, "config.csv")))

# NOTE Uncomment for resimulating
# ###
# ### draw random simulations
# ###
# using Random
# using StableRNGs
# N_nodes = 100
# N_lines = 200

# # set seed for reproducibility
# draw_task_id_seed = 3
# rng = StableRNG(draw_task_id_seed)

# # generate 100 random integers
# rand_task_ids = rand(rng, 1:length(df_config.ArrayTaskID), 100)
# rand_initial_fails = rand(rng, 1:N_lines, 100)
# seeds_rand_simulations = Dict(
#     :rand_task_ids => rand_task_ids,
#     :rand_initial_fails => rand_initial_fails)
# CSV.write(joinpath(exp_data_dir, "random_task_ids", "seeds_rand_simulations.csv"), seeds_rand_simulations)

# ###
# ### rerun simulations
# ###
# for (task_id, initial_fail) in zip(rand_task_ids, rand_initial_fails)
#     println("task_id=$task_id, initial_fail=$initial_fail")
#     sol, nw = simulate(exp_data_dir, task_id, initial_fail;
#         solverargs = (;dtmax=0.01),
#         verbose = true);
    
#     Serialization.serialize(joinpath(exp_data_dir, "random_task_ids", "sims", "task_id=$task_id,initial_fail=$initial_fail.sol"), sol)
#     Serialization.serialize(joinpath(exp_data_dir, "random_task_ids", "sims", "task_id=$task_id,initial_fail=$initial_fail.nw"), nw)
# end

###
### plot solution trajectories
###
using Colors
using CairoMakie # for normal plots
CairoMakie.activate!()


# task_id = 6048
# initial_fail = 166


seeds_rand_simulations = DataFrame(CSV.File(joinpath(exp_data_dir, "random_task_ids", "seeds_rand_simulations.csv")))
rand_task_ids = seeds_rand_simulations.rand_task_ids
rand_initial_fails = seeds_rand_simulations.rand_initial_fails

simulation_dir = joinpath(exp_data_dir, "random_task_ids", "sims")
plot_dir = joinpath(exp_data_dir, "random_task_ids", "plots")

for (task_id, initial_fail) in zip(rand_task_ids, rand_initial_fails)
    println("task_id=$task_id, initial_fail=$initial_fail")
    plot_simulation(simulation_dir, plot_dir, df_config, task_id, initial_fail)
end


# ###
# ### inspect solutions
# ###
# using NetworkDynamicsInspector
# using WGLMakie # for inspector
# # # using Bonito # for using plot pane and memorizing plots
# # # Bonito.set_cleanup_time!(720)

# task_id = 1617
# initial_fail = 15
# sol = Serialization.deserialize(joinpath(exp_data_dir, "trajectories_braess", "task_id=$task_id,initial_fail=$initial_fail.sol"));
# inspect(sol; restart=true, reset=true)
# set_sol!(sol) # optional if after inspect(sol)
# set_graphplot!(; nstate=[:ω], estate=[:S], nstate_rel=false, estate_rel=false)
# vindices = [i for i in 1:N_nodes if sol(sol.t[end], idxs=vidxs(i, :node_swing_stat))[1] != 1]
# eindices = [i for i in 1:N_lines if sol(sol.t[end], idxs=eidxs(i, :line_stat))[1] == 0]
# define_timeseries!([
#     (; selcomp=[VIndex(i) for i in vindices], states=[:ω, :ωmax], rel=false),
#     (; selcomp=[EIndex(i) for i in eindices], states=[:S, :rating], rel=false),
# ])