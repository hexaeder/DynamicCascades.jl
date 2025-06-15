#=
# Exemplary dynamic cascade in RTS-GMLC grid
![animation](rtsgmlc.mp4)
=#
# load used packages

using LinearAlgebra
print("Number of threads before setting"); print(LinearAlgebra.BLAS.get_num_threads()); print("\n")
BLAS.set_num_threads(1)
print("Number of threads after setting"); print(LinearAlgebra.BLAS.get_num_threads()); print("\n")


using Revise
using DynamicCascades
using Graphs
using MetaGraphs
using GraphMakie
using Colors
using NetworkDynamics
using Unitful
using Unitful: @u_str
using CairoMakie
using Statistics
using DataFrames
using CSV

# using GLMakie#jl
# GLMakie.activate!()#jl

damping = 0.1u"s"
scale_inertia = 1.0
network = import_system(:rtsgmlc; damping, scale_inertia, tconst = 0.01u"s")

sol = simulate(network;
        x_static=steadystate(network; tol=1e-5),
        initial_fail = Int[23],
        failtime=0.1,
        tspan = (0, 0.4),
        terminate_steady_state=true,
        trip_lines = :none,
        trip_nodes = :dynamic,
        trip_load_nodes = :none,
        f_min = -1.0/(2π),
        f_max = 1.0/(2π),
        solverargs = (;),#(;reltol=1e-7, abstol=1e-7), 
        verbose = true);

# plot solution
function plotnetwork(fig, sol, t)
    ax = Axis(fig)
    hidedecorations!(ax), hidespines!(ax)
    xlims!(ax, -10, 7)
    ylims!(ax, -5, 7)
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

fig = Figure(resolution=(1800,1000))

# plot network
fig[1,1], p = plotnetwork(fig, sol, 0.0)

# # plot all power flows
fig[2,1] = ax = Axis(fig; xlabel="time t in s", ylabel="apparent power flow in p.u.", title="power flow of all lines")
for i in 1:length(sol.load_S.saveval[1])
    t = sol.load_S.t
    y = [sol.load_S.saveval[t][i] for t in 1:length(sol.load_S.t)]
    lines!(ax, t, y; label="power flow on line ($i)", linewidth=3)
end

# plot failed power flows
# fig[2,1] = ax = Axis(fig; xlabel="time t in s", ylabel="apparent power flow in p.u.")
# for (i, l) in pairs(sol.failures.saveval)
#     t, s = seriesforidx(sol.load_S, l)
#     lines!(ax, t, s; label="flow on edge ($i)", linewidth=3)
#     scatter!(ax, (t[end], s[end]); marker=:star5, markersize=25)
# end

# plot frequencies of all gen nodes
fig[1,2] = ax = Axis(fig; xlabel="time t in s", ylabel="angular frequency ω", title="frequency transients of all generator nodes")
(nd, p, overload_cb) = nd_model(network)
state_idx = idx_containing(nd, "ω") # array: indices of ω-states
for i in state_idx
    t = sol.sol.t
    y = [sol.sol.u[t][i] for t in 1:length(sol.sol.t)]
    lines!(ax, t, y; label="frequency on node ($i)", linewidth=3)
end

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

# plot frequencies of all load nodes
fig[1,3] = ax = Axis(fig; xlabel="time t in s", ylabel="angular frequency ω", title="frequency transients of all load nodes")
for i in 1:length(sol.frequencies_load_nodes.saveval[1])
    t = sol.frequencies_load_nodes.t
    y = [sol.frequencies_load_nodes.saveval[t][i] for t in 1:length(sol.frequencies_load_nodes.t)]
    lines!(ax, t, y; label="frequency on load node ($i)", linewidth=3)
end

# plot frequencies of failed load nodes
load_node_idxs = findall(x -> x==:load, get_prop(network, 1:nv(network), :type))
fig[2,3] = ax = Axis(fig; xlabel="time t in s", ylabel="angular frequency ω", title="frequency transients of failing load nodes")
for (i, l) in pairs([findfirst(x -> x == i, load_node_idxs) for i in sol.failures_load_nodes.saveval])
    t, s = seriesforidx(sol.frequencies_load_nodes, l)
    lines!(ax, t, s; label="frequency on load node ($i)", linewidth=3)
    scatter!(ax, (t[end], s[end]); marker=:star5, markersize=25)
end

# save("250529_rtsgmlc_scaleM=$scale_inertia.pdf", fig)
fig

# Run simulation for trips on lines/nodes [23]
# Shutdown line 23 at t = 0.1
# Shutdown node 13 at t = 0.14691722754580414
# Shutdown node 23 at t = 0.3109996165841005
# Shutdown node 16 at t = 0.3183369822666648
# Shutdown node 42 at t = 0.3294000270087252
# Shutdown node 40 at t = 0.3331782029711102
# Shutdown node 45 at t = 0.33965529586705095
# Shutdown node 15 at t = 0.3399187970898225
# Shutdown node 14 at t = 0.3408124225576714
# Shutdown node 18 at t = 0.34597288817520705
# Shutdown node 21 at t = 0.34823237268629154
# Shutdown node 39 at t = 0.3506511451594256
# Shutdown node 38 at t = 0.35657436196087305
# Shutdown node 46 at t = 0.3598063781400924
# Shutdown node 7 at t = 0.36466628497929565
# Shutdown node 1 at t = 0.3652270180401851
# Shutdown node 22 at t = 0.36562566804315644
# Shutdown node 2 at t = 0.3664461672190899
# Shutdown node 71 at t = 0.3689849824257324
# Shutdown node 47 at t = 0.3707489904170369
# Shutdown node 37 at t = 0.37254016004481333
# Shutdown node 64 at t = 0.373968161941794
# Shutdown node 61 at t = 0.3743720730499023
# Shutdown node 66 at t = 0.37631385715884114
# Shutdown node 62 at t = 0.38003626784226435
# Shutdown node 69 at t = 0.38064785940096146
# Shutdown node 31 at t = 0.38087441620537005
# Shutdown node 25 at t = 0.38120893496787794
# Shutdown node 26 at t = 0.38298535762520625
# Shutdown node 63 at t = 0.3860348113011744
# Shutdown node 70 at t = 0.3919716901001555
# Shutdown node 55 at t = 0.39296970306348167
# Shutdown node 49 at t = 0.39661692220896577
# Shutdown node 50 at t = 0.3972814920117586
# ┌ Warning: Did not reach steady state! (lines none, nodes dynamic, load nodes none)