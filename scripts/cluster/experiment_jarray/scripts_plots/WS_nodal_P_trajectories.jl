"""
Watts-Strogatz-Network: Related to the curve in WS_lines+nodes_uebergang.jl with
frequency bound f_b=0.03. Plotting all frequency and power flow trajectories of all
nodes and lines that fail for I=[0.2, 3.0, 30.0].
"""

include(abspath(@__DIR__, "..", "helpers_jarray.jl"))

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
# using GLMakie
# GLMakie.activate!()

# time_stepsize = 10
#
# sol.load_S.t[1:time_stepsize:end_time]
# sol.load_S.t[1:time_stepsize:end]


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
#                tspan = (0, 100000),
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

################################################################################
#################################### debug #####################################
################################################################################
node_idx = 1
tstep = 5
sol.load_P.t[tstep]
sol.load_P.saveval[tstep] # [-0.8147757369906723, -0.8326689781607188, 0.8133352607713757, ...]
all_neighbors(network, node_idx) # 2, 77, 100
P_i = get_prop(network, node_idx, :P) # -0.06805549534093086
all_edges = collect(edges(network.graph)) # First three elements: [Edge 1 => 2, Edge 1 => 77, Edge 1 => 100, ...]
P_i + (sol.load_P.saveval[tstep][1] + sol.load_P.saveval[tstep][2] + sol.load_P.saveval[tstep][3]) == 0 # false




################################################################################
################################################################################
################################################################################

################################################################################
########################### nodal P trajectories ###############################
################################################################################

function get_ΔP_ΔS_at_node(node_idx, sol, network)
    P_i = get_prop(network, node_idx, :P)
    all_edges = collect(edges(network.graph))
    s_ΔP = Float32[] # ΔP for all time steps
    s_ΔS = Float32[] # ΔS for all time steps
    # loop over all timesteps
    for tstep in 1:length(sol.load_P.saveval)
        load_P = sol.load_P.saveval[tstep]
        load_S = sol.load_S.saveval[tstep]
        P_e = 0
        S_e = 0
        # sum over all incoming and outgoing flows at node_idx
        for j in 1:ne(network)
            #= edges are listed source to destination, with destination > source,
            e.g. for edge Edge 1 => 3 at node_idx=3 one has to multiply power flow
            by -1. =#
            if all_edges[j].src == node_idx
                P_e += load_P[j]
                #= get sign of power injection at node. this is only possible via the
                active power flows as the function `calculate_apparent_power!` returns
                only positive values =#
                S_e += (load_S[j] * sign(load_P[j]))
            end
            if all_edges[j].dst == node_idx
                P_e -= load_P[j]
                S_e -= (load_S[j] * sign(load_P[j]))
            end
        end
        # add power injection
        ΔP = P_i + P_e # TODO check
        append!(s_ΔP, ΔP)
        ΔS = P_i + S_e # TODO check
        append!(s_ΔS, ΔS)
    end
    s_ΔP, s_ΔS
end


# Solution objects
sol415 = Serialization.deserialize(joinpath(exp_data_dir, "trajectories", "task_id=415.sol"))
sol418 = Serialization.deserialize(joinpath(exp_data_dir, "trajectories", "task_id=418.sol"))
sol423 = Serialization.deserialize(joinpath(exp_data_dir, "trajectories", "task_id=423.sol"))

# 78: initial trigger line
lines_task_id_415 = [78] # sol415.failures.saveval
lines_task_id_418 = [78] # sol418.failures.saveval
lines_task_id_423 = [78, 199, 194, 131, 147, 146, 14, 106, 135, 149, 148] # sol423.failures.saveval
nodes_task_id_415 = [29, 98, 92, 58] # sol415.failures_nodes.saveval
nodes_task_id_418 = [29, 98] # sol418.failures_nodes.saveval
nodes_task_id_423 = [98, 59, 62, 61, 60] # sol423.failures_nodes.saveval

# all_failing_lines_idxs = [14, 78, 106, 131, 135, 146, 147, 148, 149, 194, 199]
# all_failing_nodes_idxs = [29, 58, 59, 60, 61, 62, 92, 98]
all_failing_lines_idxs = Int64[1:200...]
all_failing_nodes_idxs = Int64[1:100...]

node_colors = distinguishable_colors(length(all_failing_nodes_idxs)+1)
deleteat!(node_colors, 2) # delete yellow
line_colors = distinguishable_colors(length(all_failing_lines_idxs)+1)
deleteat!(line_colors, 2) # delete yellow

################################################################################
############################ Line and nodes ####################################
################################################################################
fontsize = 35
titlesize = (fontsize+5)
linewidth = 3.5
# fig = Figure(resolution=(2100,1500), fontsize= fontsize)
# xlim_30 = 36.0
fig = Figure(resolution=(3100,1500), fontsize= fontsize)
xlim_30 = .5

# FREQUENCIES ########################################################################
# frequencies of failed gen nodes I=0.2 s^2
sol = sol415
fig[1,1] = ax = Axis(fig; ylabel="Frequency [Hz]", title=L"Inertia $I=0.2$ $s^2$", titlealign = :left, titlesize = titlesize)
network = import_system_wrapper(df_config, 1)
(nd, p, overload_cb) = nd_model(network)
state_idx = idx_containing(nd, "ω") # array: indices of ω-states
node_idx = map(s -> parse(Int, String(s)[4:end]), nd.syms[state_idx]) # array: indices of vertices that are generators
# failing_nodes_idxs = sol.failures_nodes.saveval
failing_nodes_idxs = all_failing_nodes_idxs
for (i, l) in pairs([state_idx[findfirst(x -> x == i, node_idx)] for i in failing_nodes_idxs])
    t, s = seriesforidx(sol.sol, l)
    s = s./(2*π)
    node_idx = failing_nodes_idxs[i]
    lines!(ax, t, s; label="Node $node_idx", color=node_colors[i], linewidth=linewidth)
    scatter!(ax, (t[end], s[end]); color=node_colors[i], marker=:star5, markersize=25)
end
xlims!(ax, 0, 1.1)
ax.xticks = [0.0, 0.25, 0.5, 0.75, 1.0]
ax.yticks = [-0.03, -0.02, -0.01, 0.00, 0.01, 0.02, 0.03]
# axislegend(ax, position = :rt, nbanks = 3)

# frequencies of failed gen nodes I=3.0 s^2
sol = sol418
fig[1,2] = ax = Axis(fig; ylabel="Frequency [Hz]", title=L"Inertia $I=3.0$ $s^2$", titlealign = :left, titlesize = titlesize)
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
xlims!(ax, 0, 4.1)
ax.xticks = [0.1, 1, 2, 3, 4]
ax.yticks = [-0.03, -0.02, -0.01, 0.00, 0.01, 0.02, 0.03]
# axislegend(ax, position = :rt, nbanks = 2)

# frequencies of failed gen nodes I=30.0 s^2
sol = sol423
fig[1,3] = ax = Axis(fig; ylabel="Frequency [Hz]", title=L"Inertia $I=30.0$ $s^2$", titlealign = :left, titlesize = titlesize)
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
xlims!(ax, -1., xlim_30)
# ax.xticks = [0.0, 5., 10., 15., 20., 25., 30., 35.]
ax.yticks = [-0.03, -0.02, -0.01, 0.00, 0.01, 0.02, 0.03]
# axislegend(ax, position = :rt, nbanks = 2)

# FLOWS ########################################################################
# failed power flows I=0.2 s^2
sol = sol415
fig[3,1] = ax = Axis(fig; xlabel="Time [s]", ylabel="App. power flow [p.u.]")
# failing_lines_idxs = sol.failures.saveval
failing_lines_idxs = all_failing_lines_idxs
for (i, l) in pairs(failing_lines_idxs)
    t, s = seriesforidx(sol.load_S, l)
    line_idx = failing_lines_idxs[i]
    lines!(ax, t, s; label="Line $line_idx", color=line_colors[i], linewidth=linewidth)
    scatter!(ax, (t[end], s[end]); color=line_colors[i], marker=:star5, markersize=25)
end
xlims!(ax, 0, 1.1)
ax.xticks = [0.0, 0.25, 0.5, 0.75, 1.0]
# hlines!(ax, 2.1; color=:red, linestyle=:dash, linewidth=linewidth, label="Rating")
# axislegend(ax, position = (0.9,0.92), nbanks = 6, labelsize=22)

# failed power flows I=3.0 s^2
sol = sol418
fig[3,2] = ax = Axis(fig; xlabel="Time [s]", ylabel="App. power flow [p.u.]")
# failing_lines_idxs = sol.failures.saveval
failing_lines_idxs = all_failing_lines_idxs
for (i, l) in pairs(failing_lines_idxs)
    t, s = seriesforidx(sol.load_S, l)
    line_idx = failing_lines_idxs[i]
    lines!(ax, t, s; label="Line $line_idx", color=line_colors[i], linewidth=linewidth)
    scatter!(ax, (t[end], s[end]); color=line_colors[i], marker=:star5, markersize=25)
end
xlims!(ax, 0, 4.1)
# ax.xticks = [0.0, 0.1, 0.3, 0.5, 0.7, 0.9, 1.1, 1.3, 1.5]

# axislegend(ax, position = :rt, nbanks = 3)

# failed power flows I=30.0 s^2
sol = sol423
fig[3,3] = ax = Axis(fig; xlabel="Time [s]", ylabel="App. power flow [p.u.]")
# failing_lines_idxs = sol.failures.saveval
failing_lines_idxs = all_failing_lines_idxs
for (i, l) in pairs(failing_lines_idxs)
    t, s = seriesforidx(sol.load_S, l)
    line_idx = failing_lines_idxs[i]
    lines!(ax, t, s; label="Line $line_idx", color=line_colors[i], linewidth=linewidth)
    scatter!(ax, (t[end], s[end]); color=line_colors[i], marker=:star5, markersize=25)
end
xlims!(ax, -1., xlim_30)
# ax.xticks = [0.0, 5., 10., 15., 20., 25., 30., 35.]
# axislegend(ax, position = :rt, nbanks = 3)


# ΔP TRAJECTORIES ##############################################################

# ΔP trajectories I=0.2 s^2
sol = sol415
fig[2,1] = ax = Axis(fig; xlabel="Time [s]", ylabel="ΔP [p.u.]")
# failing_nodes_idxs = sol.failures_nodes.saveval
failing_nodes_idxs = all_failing_nodes_idxs
for (i, l) in pairs(failing_nodes_idxs)
    t = sol.load_P.t
    s_ΔP, s_ΔS = get_ΔP_ΔS_at_node(l, sol, network)
    node_idx = failing_nodes_idxs[i]
    lines!(ax, t, s_ΔP; label="Node $node_idx", color=node_colors[i], linewidth=linewidth)
    # plot failure time
    # tt, s = seriesforidx(sol.load_S, l)
    # scatter!(ax, (tt[end], s_ΔS[length(tt)]); color=node_colors[i], marker=:star5, markersize=25)
end
xlims!(ax, 0, 1.1)
ax.xticks = [0.0, 0.25, 0.5, 0.75, 1.0]
# hlines!(ax, 2.1; color=:red, linestyle=:dash, linewidth=linewidth, label="Rating")
# axislegend(ax, position = (0.9,0.92), nbanks = 6, labelsize=22)

# ΔP trajectories I=3.0 s^2
sol = sol418
fig[2,2] = ax = Axis(fig; xlabel="Time [s]", ylabel="ΔP [p.u.]")
# failing_nodes_idxs = sol.failures_nodes.saveval
failing_nodes_idxs = all_failing_nodes_idxs
for (i, l) in pairs(failing_nodes_idxs)
    t = sol.load_P.t
    s_ΔP, s_ΔS = get_ΔP_ΔS_at_node(l, sol, network)
    node_idx = failing_nodes_idxs[i]
    lines!(ax, t, s_ΔP; label="Node $node_idx", color=node_colors[i], linewidth=linewidth)
    # plot failure time
    # tt, s = seriesforidx(sol.load_S, l)
    # scatter!(ax, (tt[end], s_ΔS[length(tt)]); color=node_colors[i], marker=:star5, markersize=25)
end
xlims!(ax, 0, 4.1)
# ax.xticks = [0.0, 0.1, 0.3, 0.5, 0.7, 0.9, 1.1, 1.3, 1.5]

# ΔP trajectories I=30.0 s^2
sol = sol423
fig[2,3] = ax = Axis(fig; xlabel="Time [s]", ylabel="ΔP [p.u.]")
# failing_nodes_idxs = sol.failures_nodes.saveval
failing_nodes_idxs = all_failing_nodes_idxs
for (i, l) in pairs(failing_nodes_idxs)
    t = sol.load_P.t
    s_ΔP, s_ΔS = get_ΔP_ΔS_at_node(l, sol, network)
    node_idx = failing_nodes_idxs[i]
    lines!(ax, t, s_ΔS; label="Node $node_idx", color=node_colors[i], linewidth=linewidth)
    # plot failure time
    # tt, s = seriesforidx(sol.load_S, l)
    # scatter!(ax, (tt[end], s_ΔS[length(tt)]); color=node_colors[i], marker=:star5, markersize=25)
end
xlims!(ax, -1., xlim_30)
# ax.xticks = [0.0, 5., 10., 15., 20., 25., 30., 35.]

fig
# CairoMakie.save(joinpath("/home/brandner/nb_data/repos/Private_MA/paper/", "WS_traj,lines+nodes,task_ids=$task_id_array,init_fail=$initial_fail,_DeltaP_all_inertia_all_lines+nodes_long2.pdf"),fig)
# CairoMakie.save(joinpath("/home/brandner/nb_data/repos/Private_MA/paper/", "WS_traj,lines+nodes,task_ids=$task_id_array,init_fail=$initial_fail,_DeltaP_all_inertia_all_lines+nodes_long2.png"),fig)
