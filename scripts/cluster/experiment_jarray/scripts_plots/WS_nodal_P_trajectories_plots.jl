using Revise
using DynamicCascades
using Graphs
using MetaGraphs
using Unitful
using Statistics
using GraphMakie
using Colors
using DynamicCascades: PLOT_DIR

using CairoMakie

initial_fail = 78
task_id_array = [415, 418, 423]
task_id = 423
exp_name_date = "WS_k=4_exp02_PIK_HPC_K_=3,N_G=32_20240208_000237.814"

exp_data_dir = joinpath(RESULTS_DIR, exp_name_date)
df_config = DataFrame(CSV.File(joinpath(exp_data_dir, "config.csv")))
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

if steadystate_choice == :rootfind
    x_static = steadystate(network; verbose=true) # "Old" way: leads to some errors, thus the `catch`-option below
elseif steadystate_choice == :relaxation
    x_static = steadystate_relaxation(network; verbose=true) # "New" way, steady state more precise, less/no errors, probabyl slower
end

# sol = simulate(network;
#                x_static=x_static,
#                initial_fail = [initial_fail],
#                init_pert = init_pert,
#                tspan = (0, 10),
#                trip_lines = trip_lines,
#                trip_nodes = trip_nodes,
#                trip_load_nodes = :none,
#                monitored_power_flow = monitored_power_flow,
#                f_min = -freq_bound,
#                f_max = freq_bound,
#                solverargs = (;dtmax=0.01),
#                verbose = true);
#
# Serialization.serialize(joinpath(exp_data_dir, "trajectories", "task_id=$task_id.sol"), sol)


# Solution objects
sol415 = Serialization.deserialize(joinpath(exp_data_dir, "trajectories", "task_id=415.sol"))
sol418 = Serialization.deserialize(joinpath(exp_data_dir, "trajectories", "task_id=418.sol"))
sol423 = Serialization.deserialize(joinpath(exp_data_dir, "trajectories", "task_id=423.sol"))

sol = sol423

# plot solution
function plotnetwork(fig, sol, t)
    ax = Axis(fig)
    hidedecorations!(ax), hidespines!(ax)
    # xlims!(ax, -10, 7)
    # ylims!(ax, -5, 7)
    gpargs = gparguments(sol, t;
                         colortype=:abssteady,
                         Δω=Observable(2.5),
                         offlinewidth=3,
                         offlinecolor=colorant"lightgray",
                         ecolorscaling = Observable(1.0),
                         node_size=15,
                         show_labels=false)
    p = graphplot!(ax, network; gpargs...)
    return ax, p
end

# tobs = Observable(sol.sol.t[end-13])
# tobs = Observable(sol.sol.t[end])
tobs = Observable(0.0)
fig = Figure(resolution=(1800,1000))

# plot network
# fig[1,1] = Label(fig, @lift("t = "*repr(round($tobs,digits=2))*" s"), tellwidth=false)
fig[1,1], p = plotnetwork(fig, sol, tobs)

# # plot all power flows
# fig[2,1] = ax = Axis(fig; xlabel="time t in s", ylabel="apparent power flow in p.u.", title="power flow of all lines")
# for i in 1:length(sol.load_S.saveval[1])
#     t = sol.load_S.t
#     # t = sol.frequencies_load_nodes.t[1:20]
#     y = [sol.load_S.saveval[t][i] for t in 1:length(sol.load_S.t)]
#     # y = [sol.frequencies_load_nodes.saveval[t][i] for t in 1:20]
#     # lines!(ax, t, y; label="frequency on load node ($i)", linewidth=3)
#     scatter!(ax, t, y; label="power flow on line ($i)", markersize=5)
# end
# vlines!(ax, tobs; color=:black, linewidth=1)
# fig

# plot failed power flows
fig[2,1] = ax = Axis(fig; xlabel="time t in s", ylabel="apparent power flow in p.u.", title=@lift("t = "*repr(round($tobs,digits=2))*" s                                      flow transients of failing lines"))
for (i, l) in pairs(sol.failures.saveval)
    t, s = seriesforidx(sol.load_S, l)
    lines!(ax, t, s; label="flow on edge ($i)", linewidth=3)
    # scatter!(ax, t, s; label="flow on edge ($i)", markersize=5)
    scatter!(ax, (t[end], s[end]); marker=:star5, markersize=25)
end
vlines!(ax, tobs; color=:black, linewidth=1)

# plot frequencies of all gen nodes
fig[1,2] = ax = Axis(fig; xlabel="time t in s", ylabel="angular frequency ω", title="frequency transients of all generator nodes")
(nd, p, overload_cb) = nd_model(network)
state_idx = idx_containing(nd, "ω") # array: indices of ω-states
for i in state_idx
    t = sol.sol.t
    # t = sol.sol.t[1:300]
    y = [sol.sol.u[t][i] for t in 1:length(sol.sol.t)]
    # y = [sol.sol.u[t][i] for t in 1:300]
    lines!(ax, t, y; label="frequency on node ($i)", linewidth=3)
    # scatter!(ax, t, y; label="frequency on node ($i)", linewidth=3)
end
vlines!(ax, tobs; color=:black, linewidth=1)

# plot frequencies of failed gen nodes
fig[2,2] = ax = Axis(fig; xlabel="time t in s", ylabel="angular frequency ω", title="frequency transients of failing generator nodes")
(nd, p, overload_cb) = nd_model(network)
state_idx = idx_containing(nd, "ω") # array: indices of ω-states
node_idx = map(s -> parse(Int, String(s)[4:end]), nd.syms[state_idx]) # array: indices of vertices that are generators
for (i, l) in pairs([state_idx[findfirst(x -> x == i, node_idx)] for i in sol.failures_nodes.saveval])
    t, s = seriesforidx(sol.sol, l)
    lines!(ax, t, s; label="frequency on node ($i)", linewidth=3)
    scatter!(ax, (t[end], s[end]); marker=:star5, markersize=25)
end
vlines!(ax, tobs; color=:black, linewidth=1)

# plot frequencies of all load nodes
fig[1,3] = ax = Axis(fig; xlabel="time t in s", ylabel="angular frequency ω", title="frequency transients of all load nodes")
for i in 1:length(sol.frequencies_load_nodes.saveval[1])
    t = sol.frequencies_load_nodes.t
    # t = sol.frequencies_load_nodes.t[1:20]
    y = [sol.frequencies_load_nodes.saveval[t][i] for t in 1:length(sol.frequencies_load_nodes.t)]
    # y = [sol.frequencies_load_nodes.saveval[t][i] for t in 1:20]
    lines!(ax, t, y; label="frequency on load node ($i)", linewidth=3)
    # scatter!(ax, t, y; label="frequency on load node ($i)", markersize=5)
end
vlines!(ax, tobs; color=:black, linewidth=1)

# plot frequencies of failed load nodes
load_node_idxs = findall(x -> x==:load, get_prop(network, 1:nv(network), :type))
fig[2,3] = ax = Axis(fig; xlabel="time t in s", ylabel="angular frequency ω", title="frequency transients of failing load nodes")
for (i, l) in pairs([findfirst(x -> x == i, load_node_idxs) for i in sol.failures_load_nodes.saveval])
    t, s = seriesforidx(sol.frequencies_load_nodes, l)
    # scatter!(ax, t, s; label="frequency on load node ($i)", linewidth=3)
    lines!(ax, t, s; label="frequency on load node ($i)", linewidth=3)
    scatter!(ax, (t[end], s[end]); marker=:star5, markersize=25)
end
vlines!(ax, tobs; color=:black, linewidth=1)
fig


# create video
# tobs = Observable(0.0)
T = 10 #10
tmax = sol.sol.t[end] #0.15 # tmax = 3.5
# tmax = 2
tmin = 0.0
fps = 20 # 20,100
trange = range(tmin, tmax, length=Int(T * fps))

record(fig, joinpath(PLOT_DIR,"WS_2_testzz.mp4"), trange; framerate=30) do time
    tobs[] = time
end


# # single gen
# node_idx = 3
# i = state_idx[node_idx]
# t = sol.sol.t
# y = [sol.sol.u[t][i] for t in 1:length(sol.sol.t)]
# lines!(ax, t, y; label="frequency on node ($i)", linewidth=3)

# # go through all edges
# for i in 1:ne(network)
#     print("failed edge "); print(i); print("\n")
#     sol = simulate(network;
#                    initial_fail = Int[i],
#                    # tspan = (0, 50),
#                    # tspan = (0, 0.17),
#                    tspan = (0, 50),
#                    # terminate_steady_state=true,
#                    trip_lines = :dynamic,
#                    trip_nodes = :dynamic,
#                    solverargs = (;dtmax=0.01), verbose = true);
#     print("number of node failures "); print(length(sol.failures_nodes.saveval)); print("\n")
# end
