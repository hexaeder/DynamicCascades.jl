#=
# Exemplary dynamic cascade in RTS-GMLC grid
![animation](rtsgmlc.mp4)
=#
# load used packages

using Revise
# using Infiltrator
using DynamicCascades
using Graphs
using MetaGraphs
using Statistics
using GraphMakie
using Colors
using NetworkDynamics
using Unitful
using Unitful: @u_str
using DataFrames
using CSV
using CairoMakie

# using GLMakie#jl
# GLMakie.activate!()#jl

damping = 0.1u"s"
scale_inertia = 1.0
network = import_system(:rtsgmlc; damping, scale_inertia, tconst = 0.01u"s")

sol = simulate(network;
               # initial_fail = Int[37,73],
               # initial_fail = Int[27],
               initial_fail = Int[22],
               # failtime=0.01,
               tspan = (0, 50),
               # tspan = (0, 0.17),
               # tspan = (0, 0.012),
               # terminate_steady_state=true,
               trip_lines = :dynamic,
               trip_nodes = :dynamic,
               trip_load_nodes = :dynamic,
               solverargs = (;dtmax=0.01), verbose = true);

# sol.sol.u
# sol.initial_fail
# sol.failtime
# sol.trip_lines
# sol.network
# sol.load_P
# sol.failures

# plot solution
function plotnetwork(fig, sol, t)
    ax = Axis(fig)
    hidedecorations!(ax), hidespines!(ax)
    xlims!(ax, -10, 7)
    ylims!(ax, -5, 7)
    gpargs = gparguments(sol, t;
                         colortype=:abssteady,
                         Δω=Observable(2.5),
                         offlinewidth=5,
                         offlinecolor=colorant"lightgray",
                         ecolorscaling = Observable(1.0),
                         node_size=15,
                         show_labels=false)
    p = graphplot!(ax, network; gpargs...)
    return ax, p
end

tobs = Observable(sol.sol.t[end-10])
# tobs = Observable(0.0)
fig = Figure(resolution=(1800,700))

# plot network
fig[1,1] = Label(fig, @lift("t = "*repr(round($tobs,digits=2))*" s"), tellwidth=false)
fig[2,1], p = plotnetwork(fig, sol, tobs)
fig

# plot power flows
fig[2,2] = ax = Axis(fig; xlabel="time t in s", ylabel="apparent power flow in p.u.", title="flow transients of failing lines")
for (i, l) in pairs(sol.failures.saveval)
    t, s = seriesforidx(sol.load_S, l)
    lines!(ax, t, s; label="flow on edge ($i)", linewidth=3)
end
for (i, l) in pairs(sol.failures.saveval)
    t, s = seriesforidx(sol.load_S, l)
    scatter!(ax, (t[end], s[end]); marker=:star5, markersize=25)
end
vlines!(ax, tobs; color=:black, linewidth=1)

# # plot frequencies of failed gen nodes
# fig[2,3] = ax = Axis(fig; xlabel="time t in s", ylabel="frequency", title="frequency transients of failing gen nodes")
# (nd, p, overload_cb) = nd_model(network)
# state_idx = idx_containing(nd, "ω") # array: indices of ω-states
# node_idx = map(s -> parse(Int, String(s)[4:end]), nd.syms[state_idx]) # array: indices of vertices that are generators
#
# state_idx_failures = []
# for i in sol.failures_nodes.saveval
#     push!(state_idx_failures, state_idx[findfirst(x -> x == i, node_idx)])
# end
# for (i, l) in pairs(state_idx_failures)
#     t, s = seriesforidx(sol.sol, l)
#     lines!(ax, t, s; label="frequency on node ($i)", linewidth=3)
# end
# for (i, l) in pairs(state_idx_failures)
#     t, s = seriesforidx(sol.sol, l)
#     scatter!(ax, (t[end], s[end]); marker=:star5, markersize=25)
#     # print(s[end]); print("\n")
# end
# vlines!(ax, tobs; color=:black, linewidth=1)
# fig

# # plot frequencies of...
# # all gen nodes
# fig[2,3] = ax = Axis(fig; xlabel="time t in s", ylabel="frequency", title="frequency transients of all gen nodes")
# (nd, p, overload_cb) = nd_model(network)
# state_idx = idx_containing(nd, "ω") # array: indices of ω-states
#
# for i in state_idx
#     t = sol.sol.t
#     # t = sol.sol.t[1:300]
#     y = [sol.sol.u[t][i] for t in 1:length(sol.sol.t)]
#     # y = [sol.sol.u[t][i] for t in 1:300]
#     lines!(ax, t, y; label="frequency on node ($i)", linewidth=5)
# end
# fig
# save(joinpath(PLOT_DIR, "rtsgmlc_scaleM=$scale_inertia.pdf"), fig)

# plot frequencies of...
# all load nodes
fig[2,3] = ax = Axis(fig; xlabel="time t in s", ylabel="frequency", title="frequency transients of all load nodes")
for i in 1:length(sol.frequencies_load_nodes.saveval[1])
    t = sol.frequencies_load_nodes.t
    y = [sol.frequencies_load_nodes.saveval[t][i] for t in 1:length(sol.frequencies_load_nodes.t)]
    lines!(ax, t, y; label="frequency on load node ($i)", linewidth=5)
end
fig

# # single gen
# node_idx = 3
# i = state_idx[node_idx]
# t = sol.sol.t
# y = [sol.sol.u[t][i] for t in 1:length(sol.sol.t)]
# lines!(ax, t, y; label="frequency on node ($i)", linewidth=5)

T = 20
tmax = sol.sol.t[end] #0.15 # tmax = 3.5
# tmax = 3.0
tmin = 0.0
fps = 30
# fps = 100
trange = range(tmin, tmax, length=Int(T * fps))

record(fig, joinpath(PLOT_DIR,"rtsgmlc_scaleM=$scale_inertia.mp4"), trange; framerate=30) do time
    tobs[] = time
end

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
