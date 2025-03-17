"""
"""

include(abspath(@__DIR__, "..", "helpers_jarray.jl"))
include(abspath(@__DIR__, "WS_trajectories_new_ND_single_model_port.jl"))

using DynamicCascades
# using Graphs
using MetaGraphs
# using Unitful
# using Statistics
# using GraphMakie
# using Colors
# using DynamicCascades: PLOT_DIR

# using CairoMakie



###
### read in parameters from .csv
###
initial_fail = 78
task_id = 1720
task_id_array = [1720, 1723, 1728]

# exp_name_date = "WS_k=4_exp02_PIK_HPC_K_=3,N_G=32_20240208_000237.814"
exp_name_date = "WS_k=4_exp04_vary_I_only_lines_and_nodes_PIK_HPC_K_=3,N_G=32_20250126_012357.344"
exp_data_dir = joinpath(RESULTS_DIR, exp_name_date)
df_config = DataFrame(CSV.File(joinpath(exp_data_dir, "config.csv")))
N,k,β,graph_seed,μ,σ,distr_seed,K,α,M,γ,τ,freq_bound,trip_lines,trip_nodes,init_pert,ensemble_element = get_network_args_stripped(df_config, task_id)
network = import_system_wrapper(df_config, task_id)


###
### build NetworkDynamics.jl Network
###

# loop over vertices and assign vertex models & parameter
vm_array = VertexModel[]
for i in 1:nv(network)
    # P = P_inj - P_load see `balance_power!`
    # P_inj = P + P_load
    P = get_prop(network, i, :P)
    Pload = get_prop(network, i, :P_load)
    Pmech = P + Pload

    type = get_prop(network, i, :type)
    if type == :gen
        vm = SwingDynLoadModel(M=M,D=γ,τ=τ,ωmax=freq_bound,Pmech=Pmech,Pload=Pload)
    elseif type == :load
        vm = DynLoadModel(τ=τ,Pload=Pload)  
    end
    push!(vm_array, vm)
end

# generate `Network` object
nw = Network(network.graph, vm_array, Line(K=K,rating=α*K); dealias=true)

# Check if network is power balanced
nw_state = NWState(nw)
p = nw_state.p
@assert isapprox(sum(p.v[1:100, :Pload]), sum(p.v[1:100, :Pmech]))

# set initial perturbation CB
init_perturb = PresetTimeComponentCallback(0.1,
    ComponentAffect([], [:active]) do u, p, ctx
        println("Shutdown line $(ctx.eidx) at t=$(ctx.integrator.t)")
        p[:active] = 0
    end
)
set_callback!(nw[EIndex(initial_fail)], init_perturb)


s0=find_fixpoint(nw)
prob = ODEProblem(nw, uflat(s0), (0, 1), pflat(s0), callback=get_callbacks(nw)); # TODO warum `pflat(s0)` warum extrahiert man Parameter  nicht aus `nw`?
sol = solve(prob, Rodas4P())



####
network.graph
collect(edges(network.graph))

# x0=NWState(nw)
# p=x0.p
# x_static0 = solve(NonlinearProblem(nw, x0, p), NLSolveJL())

s0 = find_fixpoint(nw, x0, p) # TODO check if I get the same fixpoint.
prob = ODEProblem(nw, uflat(s0), (0, 10), pflat(s0), callback=get_callbacks(nw)); # TODO warum `pflat(s0)` warum extrahiert man Parameter  nicht aus `nw`?
sol = solve(prob, Tsit5())
####

# call for interactive inspection
# inspect(sol)

let
    fig = Figure(size=(800, 800))
    ax = Axis(fig[1, 1]; xlabel="time", ylabel="voltage angle θ")
    lines!(ax, sol, idxs=vidxs(1:2, :θ))
    axislegend(ax)
    ax = Axis(fig[2, 1]; xlabel="time", ylabel="rotor frequency ω")
    lines!(ax, sol, idxs=vidxs(2, :ω))
    axislegend(ax)
    ax = Axis(fig[3, 1]; xlabel="time", ylabel="injected power P")
    lines!(ax, sol, idxs=vidxs(1:2, :Pinj))
    axislegend(ax)
    fig
end



























#############################################################################################
#############################################################################################
#############################################################################################

#= Generate and save sol-Objects, the data for the plots. Only use this for recreating
the solution objects =#
####
#### save and load sol objects
####
# Serialization.serialize(joinpath(exp_data_dir, "trajectories", "task_id=$task_id.sol"), sol) TODO add initially failed line

# # Solution objects
# sol415 = Serialization.deserialize(joinpath(exp_data_dir, "trajectories", "task_id=415.sol"))
# sol418 = Serialization.deserialize(joinpath(exp_data_dir, "trajectories", "task_id=418.sol"))
# sol423 = Serialization.deserialize(joinpath(exp_data_dir, "trajectories", "task_id=423.sol"))



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

# CairoMakie.save(joinpath(MA_DIR, "WS", "WS_traj,lines_only,task_ids=$task_id_array,init_fail=$initial_fail.pdf"),fig)
# CairoMakie.save(joinpath(MA_DIR, "WS", "WS_traj,lines_only,task_ids=$task_id_array,init_fail=$initial_fail.png"),fig)


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
