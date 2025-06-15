using Revise

include(abspath(@__DIR__, "..", "helpers_jarray.jl"))

using DynamicCascades
using Graphs
using MetaGraphs
using Unitful
using Statistics
using GraphMakie
using Colors
using DynamicCascades: PLOT_DIR

using CairoMakie


initial_fail = 95
task_id = 2310
exp_name_date = "WS_k=4_exp04_vary_I_only_lines_and_nodes_PIK_HPC_K_=3,N_G=32_20250321_171511.976"

exp_data_dir = joinpath(RESULTS_DIR, exp_name_date)
df_config = DataFrame(CSV.File(joinpath(exp_data_dir, "config.csv")))
# adjust filepaths 
df_config[!, :filepath_graph] = replace.(df_config[!, :filepath_graph],"/home/brandner" => "/home/brandner/nb_data/HU_Master/2122WS/MA")
df_config[!, :filepath_steady_state] = replace.(df_config[!, :filepath_steady_state],"/home/brandner" => "/home/brandner/nb_data/HU_Master/2122WS/MA")
exp_params_dict = Serialization.deserialize(joinpath(exp_data_dir, "exp.params"))

N,k,β,graph_seed,μ,σ,distr_seed,K,α,M,γ,τ,freq_bound,trip_lines,trip_nodes,init_pert,ensemble_element = get_network_args_stripped(df_config, task_id)
monitored_power_flow = exp_params_dict[:monitored_power_flow]
steadystate_choice = exp_params_dict[:steadystate_choice]
# steadystate_choice = :relaxation

#= Generate and save sol-Objects, the data for the plots. Only use this for recreating
the solution objects =#
network = import_system_wrapper(df_config, task_id)

# network is balanced
sum([get_prop(network, i, :P) for i in 1:nv(network)]) # 3.3306690738754696e-15

# read in steady state
steady_state_dict  = CSV.File(df_config[task_id,:filepath_steady_state])
x_static = steady_state_dict[:SteadyState]

node_failure_model = :no_failures

if node_failure_model == :no_failures
    trip_lines = :none
    trip_nodes = :none
end

sol = simulate(network;
               x_static=x_static,
               initial_fail = [initial_fail],
               init_pert = init_pert,
               tspan = (0, 15),
               trip_lines = trip_lines,
               trip_nodes = trip_nodes,
               trip_load_nodes = :none,
               node_failure_model = node_failure_model,
               monitored_power_flow = monitored_power_flow,
               f_min = -freq_bound,
               f_max = freq_bound,
               solverargs = (;dtmax=0.01),
               verbose = true);
               
Serialization.serialize(joinpath(exp_data_dir, "testing_new_node_failure_models_for_old_ND_code", "node_failure_model=$node_failure_model,task_id=$task_id,initial_fail=$initial_fail.sol"), sol)

# calculate indices of failing lines and nodes
all_failing_nodes_idxs = sol.failures_nodes.saveval
all_failing_lines_idxs = sol.failures.saveval
node_colors_failing = distinguishable_colors(length(all_failing_nodes_idxs))
line_colors_failing = distinguishable_colors(length(all_failing_lines_idxs))

# indices of all lines and nodes
all_nodes_idxs = [i for i in 1:nv(network)]
all_lines_idxs = [i for i in 1:ne(network)]
node_colors = distinguishable_colors(length(all_nodes_idxs))
line_colors = distinguishable_colors(length(all_lines_idxs))

################################################################################
############################ Line and nodes ####################################
################################################################################
fontsize = 35
titlesize = (fontsize+5)
linewidth = 3.5
fig = Figure(resolution=(3100,1500), fontsize=fontsize)
# Add a global title in the first row spanning all columns.
xlim = sol.sol.t[end]

N,k,β,graph_seed,μ,σ,distr_seed,K,α,M,γ,τ,freq_bound,trip_lines,trip_nodes,init_pert,ensemble_element = get_network_args_stripped(df_config, task_id)

# FREQUENCIES ########################################################################
fig[1,1] = ax = Axis(fig; ylabel="Frequency [Hz]", title="Lines & nodes that fail (left column), all lines & nodes (right column):  I=$M,D=$γ,τ=$τ,f_b=$freq_bound,α=$α,K=$K,N=$N,k=$k,β=$β,ensemble element=$ensemble_element", titlealign = :left, titlesize = titlesize)

network = import_system_wrapper(df_config, 1)
(nd, p, overload_cb) = nd_model(network)
state_idx = idx_containing(nd, "ω") # array: indices of ω-states
node_idx = map(s -> parse(Int, String(s)[4:end]), nd.syms[state_idx]) # array: indices of vertices that are generators
failing_nodes_idxs = all_failing_nodes_idxs
for (i, l) in pairs([state_idx[findfirst(x -> x == i, node_idx)] for i in failing_nodes_idxs])
    t, s = seriesforidx(sol.sol, l)
    s = s./(2*π)
    node_idx = failing_nodes_idxs[i]
    lines!(ax, t, s; label="Node $node_idx", color=node_colors[i], linewidth=linewidth)
    scatter!(ax, (t[end], s[end]); color=node_colors[i], marker=:star5, markersize=25)
end

fig[1,2] = ax = Axis(fig; titlealign = :left, titlesize = titlesize)
(nd, p, overload_cb) = nd_model(network)
state_idx = idx_containing(nd, "ω") # array: indices of ω-states
node_idx = map(s -> parse(Int, String(s)[4:end]), nd.syms[state_idx]) # array: indices of vertices that are generators
failing_nodes_idxs = all_nodes_idxs
for (i, l) in pairs([state_idx[findfirst(x -> x == i, node_idx)] for i in failing_nodes_idxs])
    t, s = seriesforidx(sol.sol, l)
    s = s./(2*π)
    node_idx = failing_nodes_idxs[i]
    lines!(ax, t, s; label="Node $node_idx", color=node_colors[i], linewidth=linewidth)
    scatter!(ax, (t[end], s[end]); color=node_colors[i], marker=:star5, markersize=25)
end
xlims!(ax, 0, xlim)

# FLOWS ########################################################################
fig[2,1] = ax = Axis(fig; xlabel="Time [s]", ylabel="App. power flow [p.u.]")
failing_lines_idxs = all_failing_lines_idxs
for (i, l) in pairs(failing_lines_idxs)
    t, s = seriesforidx(sol.load_S, l)
    line_idx = failing_lines_idxs[i]
    lines!(ax, t, s; label="Line $line_idx", color=line_colors[i], linewidth=linewidth)
    scatter!(ax, (t[end], s[end]); color=line_colors[i], marker=:star5, markersize=25)
end

fig[2,2] = ax = Axis(fig; xlabel="Time [s]")
failing_lines_idxs = all_lines_idxs
for (i, l) in pairs(failing_lines_idxs)
    t, s = seriesforidx(sol.load_S, l)
    line_idx = failing_lines_idxs[i]
    lines!(ax, t, s; label="Line $line_idx", color=line_colors[i], linewidth=linewidth)
    scatter!(ax, (t[end], s[end]); color=line_colors[i], marker=:star5, markersize=25)
end
xlims!(ax, 0, xlim)
CairoMakie.save(joinpath(exp_data_dir, "testing_new_node_failure_models_for_old_ND_code", "node_failure_model=$node_failure_model,ensemble_element=$ensemble_element,I=$M,D=$γ,f_b=$freq_bound,task_id=$task_id,initial_fail=$initial_fail.pdf"),fig)
CairoMakie.save(joinpath(exp_data_dir, "testing_new_node_failure_models_for_old_ND_code", "node_failure_model=$node_failure_model,ensemble_element=$ensemble_element,I=$M,D=$γ,f_b=$freq_bound,task_id=$task_id,initial_fail=$initial_fail.png"),fig)
fig

