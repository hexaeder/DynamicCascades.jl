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
task_id = 418
exp_name_date = "WS_k=4_exp02_PIK_HPC_K_=3,N_G=32_20240208_000237.814"

exp_data_dir = joinpath(RESULTS_DIR, exp_name_date)
df_config = DataFrame(CSV.File(joinpath(exp_data_dir, "config.csv")))
exp_params_dict = Serialization.deserialize(joinpath(exp_data_dir, "exp.params"))

N,k,β,graph_seed,μ,σ,distr_seed,K,α,M,γ,τ,freq_bound,trip_lines,trip_nodes,init_pert,ensemble_element = get_network_args_stripped(df_config, task_id)
monitored_power_flow = exp_params_dict[:monitored_power_flow]
steadystate_choice = exp_params_dict[:steadystate_choice]

#= Generate and save sol-Objects, the data for the plots. Only use this for recreating
the solution objects =#
# network = import_system_wrapper(df_config, task_id)
#
# if steadystate_choice == :rootfind
#     x_static = steadystate(network; verbose=true) # "Old" way: leads to some errors, thus the `catch`-option below
# elseif steadystate_choice == :relaxation
#     x_static = steadystate_relaxation(network; verbose=true) # "New" way, steady state more precise, less/no errors, probabyl slower
# end
#
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

all_failing_lines_idxs = [14, 78, 106, 131, 135, 146, 147, 148, 149, 194, 199]
all_failing_nodes_idxs = [29, 58, 59, 60, 61, 62, 92, 98]

node_colors = distinguishable_colors(9)
deleteat!(node_colors, 2) # delete yellow
line_colors = distinguishable_colors(12)
deleteat!(line_colors, 2) # delete yellow


################################################################################
############################ Line and nodes ####################################
################################################################################
fontsize = 26
titlesize = (fontsize+5)
linewidth = 3.5
fig = Figure(resolution=(2100,1500), fontsize= fontsize)

# FREQUENCIES ########################################################################
# frequencies of failed gen nodes I=0.2 s^2
sol = sol415
fig[1,1] = ax = Axis(fig; title=L"Inertia $I=0.2$ $s^2$", titlealign = :left, titlesize = titlesize)
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
xlims!(ax, 0, 1.3)
ax.xticks = [0.0, 0.1, 0.3, 0.5, 0.7, 0.9, 1.1]
ax.yticks = [-0.03, -0.02, -0.01, 0.00, 0.01, 0.02, 0.03]
axislegend(ax, position = :rt, nbanks = 2)

# frequencies of failed gen nodes I=3.0 s^2
sol = sol418
fig[2,1] = ax = Axis(fig; ylabel="Frequency [Hz]", title=L"Inertia $I=3.0$ $s^2$", titlealign = :left, titlesize = titlesize)
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
xlims!(ax, 0, 5.0)
ax.xticks = [0.1, 1, 2, 3, 4]
ax.yticks = [-0.03, -0.02, -0.01, 0.00, 0.01, 0.02, 0.03]
# axislegend(ax, position = :rt, nbanks = 2)

# frequencies of failed gen nodes I=30.0 s^2
sol = sol423
fig[3,1] = ax = Axis(fig; xlabel="Time [s]", title=L"Inertia $I=30.0$ $s^2$", titlealign = :left, titlesize = titlesize)
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
xlims!(ax, 0, 75.0)
ax.xticks = [0, 25, 50, 75, 100]
ax.yticks = [-0.03, -0.02, -0.01, 0.00, 0.01, 0.02, 0.03]
# axislegend(ax, position = :rt, nbanks = 2)

# FLOWS ########################################################################
# failed power flows I=0.2 s^2
sol = sol415
fig[1,2] = ax = Axis(fig; title=L"Inertia $I=0.2$ $s^2$", titlealign = :left, titlesize = titlesize)
# failing_lines_idxs = sol.failures.saveval
failing_lines_idxs = all_failing_lines_idxs
for (i, l) in pairs(failing_lines_idxs)
    t, s = seriesforidx(sol.load_S, l)
    line_idx = failing_lines_idxs[i]
    lines!(ax, t, s; label="Line $line_idx", color=line_colors[i], linewidth=linewidth)
    scatter!(ax, (t[end], s[end]); color=line_colors[i], marker=:star5, markersize=25)
end
xlims!(ax, 0, 1.1)
ax.xticks = [0.0, 0.1, 0.3, 0.5, 0.7, 0.9, 1.1]
# hlines!(ax, 2.1; color=:red, linestyle=:dash, linewidth=linewidth, label="Rating")
axislegend(ax, position = (0.9,0.92), nbanks = 6, labelsize=22)

# failed power flows I=3.0 s^2
sol = sol418
fig[2,2] = ax = Axis(fig; ylabel="Active power flow [p.u.]", title=L"Inertia $I=3.0$ $s^2$", titlealign = :left, titlesize = titlesize)
# failing_lines_idxs = sol.failures.saveval
failing_lines_idxs = all_failing_lines_idxs
for (i, l) in pairs(failing_lines_idxs)
    t, s = seriesforidx(sol.load_S, l)
    line_idx = failing_lines_idxs[i]
    lines!(ax, t, s; label="Line $line_idx", color=line_colors[i], linewidth=linewidth)
    scatter!(ax, (t[end], s[end]); color=line_colors[i], marker=:star5, markersize=25)
end
xlims!(ax, 0, 1.5)
ax.xticks = [0.0, 0.1, 0.3, 0.5, 0.7, 0.9, 1.1, 1.3, 1.5]

# axislegend(ax, position = :rt, nbanks = 3)

# failed power flows I=30.0 s^2
sol = sol423
fig[3,2] = ax = Axis(fig; xlabel="Time [s]", title=L"Inertia $I=30.0$ $s^2$", titlealign = :left, titlesize = titlesize)
# failing_lines_idxs = sol.failures.saveval
failing_lines_idxs = all_failing_lines_idxs
for (i, l) in pairs(failing_lines_idxs)
    t, s = seriesforidx(sol.load_S, l)
    line_idx = failing_lines_idxs[i]
    lines!(ax, t, s; label="Line $line_idx", color=line_colors[i], linewidth=linewidth)
    scatter!(ax, (t[end], s[end]); color=line_colors[i], marker=:star5, markersize=25)
end
# xlims!(ax, 0, 1.7)
# ax.xticks = [0.0, 0.1, 0.3, 0.5, 0.7, 0.9, 1.1]
# axislegend(ax, position = :rt, nbanks = 3)
fig
CairoMakie.save(joinpath(MA_DIR, "WS", "WS_traj,lines+nodes,task_ids=$task_id_array,init_fail=$initial_fail.pdf"),fig)
CairoMakie.save(joinpath(MA_DIR, "WS", "WS_traj,lines+nodes,task_ids=$task_id_array,init_fail=$initial_fail.png"),fig)


################################################################################
############################ Nodes only ########################################
################################################################################

fontsize = 26
titlesize = (fontsize+5)
linewidth = 3.5
fig = Figure(resolution=(1000,1300), fontsize= fontsize)
# fig = Figure(fontsize= fontsize)
# FREQUENCIES ########################################################################
# frequencies of failed gen nodes I=0.2 s^2
sol = sol415
fig[1,1] = ax = Axis(fig; title=L"Inertia $I=0.2$ $s^2$", titlealign = :left, titlesize = titlesize)
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
xlims!(ax, 0, 1.3)
ax.xticks = [0.0, 0.1, 0.3, 0.5, 0.7, 0.9, 1.1]
ax.yticks = [-0.03, -0.02, -0.01, 0.00, 0.01, 0.02, 0.03]
axislegend(ax, position = :rt, nbanks = 2)

# frequencies of failed gen nodes I=3.0 s^2
sol = sol418
fig[2,1] = ax = Axis(fig; ylabel="Frequency [Hz]", title=L"Inertia $I=3.0$ $s^2$", titlealign = :left, titlesize = titlesize)
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
xlims!(ax, 0, 5.0)
ax.xticks = [0.1, 1, 2, 3, 4]
ax.yticks = [-0.03, -0.02, -0.01, 0.00, 0.01, 0.02, 0.03]
# axislegend(ax, position = :rt, nbanks = 2)

# frequencies of failed gen nodes I=30.0 s^2
sol = sol423
fig[3,1] = ax = Axis(fig; xlabel="Time [s]", title=L"Inertia $I=30.0$ $s^2$", titlealign = :left, titlesize = titlesize)
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
xlims!(ax, 0, 75.0)
ax.xticks = [0, 25, 50, 75, 100]
ax.yticks = [-0.03, -0.02, -0.01, 0.00, 0.01, 0.02, 0.03]
# axislegend(ax, position = :rt, nbanks = 2)

CairoMakie.save(joinpath(MA_DIR, "WS", "WS_traj,nodes_only,task_ids=$task_id_array,init_fail=$initial_fail.pdf"),fig)
CairoMakie.save(joinpath(MA_DIR, "WS", "WS_traj,nodes_only,task_ids=$task_id_array,init_fail=$initial_fail.png"),fig)


################################################################################
############################ Lines only ########################################
################################################################################

fontsize = 26
titlesize = (fontsize+5)
linewidth = 3.5
fig = Figure(resolution=(1100,1400), fontsize= fontsize)
# FLOWS ########################################################################
# failed power flows I=0.2 s^2
sol = sol415
fig[1,1] = ax = Axis(fig; title=L"Inertia $I=0.2$ $s^2$", titlealign = :left, titlesize = titlesize)
# failing_lines_idxs = sol.failures.saveval
failing_lines_idxs = all_failing_lines_idxs
for (i, l) in pairs(failing_lines_idxs)
    t, s = seriesforidx(sol.load_S, l)
    line_idx = failing_lines_idxs[i]
    lines!(ax, t, s; label="Line $line_idx", color=line_colors[i], linewidth=linewidth)
    scatter!(ax, (t[end], s[end]); color=line_colors[i], marker=:star5, markersize=25)
end
xlims!(ax, 0, 1.1)
ax.xticks = [0.0, 0.1, 0.3, 0.5, 0.7, 0.9, 1.1]
# hlines!(ax, 2.1; color=:red, linestyle=:dash, linewidth=linewidth, label="Rating")
axislegend(ax, position = (0.9,0.93), nbanks = 6, labelsize=22)

# failed power flows I=3.0 s^2
sol = sol418
fig[2,1] = ax = Axis(fig; ylabel="Active power flow [p.u.]", title=L"Inertia $I=3.0$ $s^2$", titlealign = :left, titlesize = titlesize)
# failing_lines_idxs = sol.failures.saveval
failing_lines_idxs = all_failing_lines_idxs
for (i, l) in pairs(failing_lines_idxs)
    t, s = seriesforidx(sol.load_S, l)
    line_idx = failing_lines_idxs[i]
    lines!(ax, t, s; label="Line $line_idx", color=line_colors[i], linewidth=linewidth)
    scatter!(ax, (t[end], s[end]); color=line_colors[i], marker=:star5, markersize=25)
end
xlims!(ax, 0, 1.5)
ax.xticks = [0.0, 0.1, 0.3, 0.5, 0.7, 0.9, 1.1, 1.3, 1.5]

# axislegend(ax, position = :rt, nbanks = 3)

# failed power flows I=30.0 s^2
sol = sol423
fig[3,1] = ax = Axis(fig; xlabel="Time [s]", title=L"Inertia $I=30.0$ $s^2$", titlealign = :left, titlesize = titlesize)
# failing_lines_idxs = sol.failures.saveval
failing_lines_idxs = all_failing_lines_idxs
for (i, l) in pairs(failing_lines_idxs)
    t, s = seriesforidx(sol.load_S, l)
    line_idx = failing_lines_idxs[i]
    lines!(ax, t, s; label="Line $line_idx", color=line_colors[i], linewidth=linewidth)
    scatter!(ax, (t[end], s[end]); color=line_colors[i], marker=:star5, markersize=25)
end
# xlims!(ax, 0, 1.7)
# ax.xticks = [0.0, 0.1, 0.3, 0.5, 0.7, 0.9, 1.1]
# axislegend(ax, position = :rt, nbanks = 3)

CairoMakie.save(joinpath(MA_DIR, "WS", "WS_traj,lines_only,task_ids=$task_id_array,init_fail=$initial_fail.pdf"),fig)
CairoMakie.save(joinpath(MA_DIR, "WS", "WS_traj,lines_only,task_ids=$task_id_array,init_fail=$initial_fail.png"),fig)


# # plot all power flows # NOTE This takes really long
# fig[2,2] = ax = Axis(fig; xlabel="time t in s", ylabel="apparent power flow in p.u.", title="power flow of all lines")
# for i in 1:length(sol.load_S.saveval[1])
#     if t_end == :plot_all_timesteps
#         t = sol.load_S.t
#     else
#         t = sol.load_S.t[1:1:t_end]
#     end
#
#     # if M > 5.0
#     #     t = sol.load_S.t[1:time_stepsize:end]
#     # end
#     # t = sol.frequencies_load_nodes.t[1:20]
#     y = [sol.load_S.saveval[j][i] for j in 1:length(t)]
#     # y = [sol.frequencies_load_nodes.saveval[t][i] for t in 1:20]
#     lines!(ax, t, y; label="frequency on load node ($i)", linewidth=2)
#     # scatter!(ax, t, y; label="power flow on line ($i)", markersize=5)
# end
# vlines!(ax, tobs; color=:black, linewidth=1)

# # plot frequencies of all gen nodes
# fig[1,2] = ax = Axis(fig; xlabel="time t in s", ylabel="angular frequency ω", title=@lift("t = "*repr(round($tobs,digits=2))*" s        frequency transients of all generator nodes"))
# (nd, p, overload_cb) = nd_model(network)
# state_idx = idx_containing(nd, "ω") # array: indices of ω-states
# for i in state_idx
#     if t_end == :plot_all_timesteps
#         t = sol.sol.t
#     else
#         t = sol.sol.t[1:1:t_end]
#     end
#     # t = sol.sol.t[1:1:t_end]
#     # if M > 5.0
#     #     t = sol.sol.t[1:time_stepsize:end]
#     # end
#     y = [sol.sol.u[j][i] for j in 1:length(t)]
#     # y = [sol.sol.u[t][i] for t in 1:300]
#     lines!(ax, t, y; label="frequency on node ($i)", linewidth=2)
#     # scatter!(ax, t, y; label="frequency on node ($i)", linewidth=3)
# end
# vlines!(ax, tobs; color=:black, linewidth=1)


# # plot frequencies of all load nodes
# fig[1,3] = ax = Axis(fig; xlabel="time t in s", ylabel="angular frequency ω", title="frequency transients of all load nodes")
# for i in 1:length(sol.frequencies_load_nodes.saveval[1])
#     t = sol.frequencies_load_nodes.t
#     # t = sol.frequencies_load_nodes.t[1:20]
#     y = [sol.frequencies_load_nodes.saveval[t][i] for t in 1:length(sol.frequencies_load_nodes.t)]
#     # y = [sol.frequencies_load_nodes.saveval[t][i] for t in 1:20]
#     lines!(ax, t, y; label="frequency on load node ($i)", linewidth=2)
#     # scatter!(ax, t, y; label="frequency on load node ($i)", markersize=5)
# end
# vlines!(ax, tobs; color=:black, linewidth=1)

# # plot frequencies of failed load nodes
# load_node_idxs = findall(x -> x==:load, get_prop(network, 1:nv(network), :type))
# fig[2,3] = ax = Axis(fig; xlabel="time t in s", ylabel="angular frequency ω", title="frequency transients of failing load nodes")
# for (i, l) in pairs([findfirst(x -> x == i, load_node_idxs) for i in sol.failures_load_nodes.saveval])
#     t, s = seriesforidx(sol.frequencies_load_nodes, l)
#     # scatter!(ax, t, s; label="frequency on load node ($i)", linewidth=3)
#     lines!(ax, t, s; label="frequency on load node ($i)", linewidth=3)
#     scatter!(ax, (t[end], s[end]); marker=:star5, markersize=25)
# end
# vlines!(ax, tobs; color=:black, linewidth=1)

# # create video
# T = 20 #10
# tmin = 0.0
# if t_end == :plot_all_timesteps
#     tmax = sol.sol.t[end]
# else
#     tmax = sol.sol.t[t_end]
# end
#
# fps = 20 # 20,100
# trange = range(tmin, tmax, length=Int(T * fps))
#
# string_args = string_network_args(df_config, task_id)
# record(fig, joinpath(MA_DIR, "WS", "WS_traj,t_end=$t_end,task_id=$task_id,init_fail=$initial_fail,$string_args.mp4"), trange; framerate=30) do time
#     tobs[] = time
# end
