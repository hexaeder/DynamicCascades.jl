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
scale_inertia = 1.1
network = import_system(:rtsgmlc; damping, scale_inertia, tconst = 0.01u"s")

# network = import_system(:rtsgmlc; damping, scale_inertia, tconst = 0.11u"s")
# network = import_system(:rtsgmlc; damping, scale_inertia, tconst = 0.5u"s")

sol = simulate(network;
               # initial_fail = Int[37,73],
               # initial_fail = Int[27],
               # initial_fail = Int[27],
               initial_fail = Int[11],
               failtime=0.025,
               init_pert = :line,
               # P_perturb = 1.0,
               # tspan = (0, 0.02),
               # tspan = (0, 0.05),
               # tspan = (0, 0.17),
               # tspan = (0, 0.012),
               # tspan = (0, 0.4),
               # tspan = (0, 0.11),
               # tspan = (0, 0.08),
               tspan = (0, 3),
               terminate_steady_state=true,
               trip_lines = :dynamic,
               trip_nodes = :dynamic,
               trip_load_nodes = :none,
               # f_min = -1.0/(2π),
               # f_max = 1.0/(2π),
               f_min = -2.5,
               f_max = 1.5,
               solverargs = (;),#(;reltol=1e-7, abstol=1e-7),
               verbose = true);
               # solverargs = (;dtmax=0.01, dtmin=0.0001, force_dtmin=true, save_everystep=false), verbose = true);
               # solverargs = (;dtmax=0.0000001), verbose = true);

main()
sol.sol.t
# save(joinpath(PLOT_DIR, "rtsgmlc_scaleM=$scale_inertia.pdf"), fig)

# sol.sol.u
# sol.initial_fail
# sol.failtime
# sol.trip_lines
# sol.network
# sol.load_P
# sol.failures

# nd, = nd_model(network)
# ω_state_idxs = idx_containing(nd, "ω")
# gen_node_idxs = map(s -> parse(Int, String(s)[4:end]), nd.syms[ω_state_idxs])

function main()

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
end

# create video
# tobs = Observable(0.0)
T = 20 #10
tmax = sol.sol.t[end] #0.15 # tmax = 3.5
# tmax = 2
tmin = 0.0
fps = 30 # 20,100
trange = range(tmin, tmax, length=Int(T * fps))

record(fig, joinpath(PLOT_DIR,"test_init_node_corrected_gen_node_fails_included_rtsgmlc_scaleM=$scale_inertia.mp4"), trange; framerate=30) do time
    tobs[] = time
end


# single gen
node_idx = 3
i = state_idx[node_idx]
t = sol.sol.t
y = [sol.sol.u[t][i] for t in 1:length(sol.sol.t)]
lines!(ax, t, y; label="frequency on node ($i)", linewidth=3)

# go through all/certain edges
edges_to_fail = [22,23,33,59]
# edges_to_fail = [22]
for i in edges_to_fail
    print("failed edge "); print(i); print("\n")
    sol = simulate(network;
                   initial_fail = Int[i],
                   init_pert = :line,
                   # tspan = (0, 50),
                   # tspan = (0, 0.17),
                   tspan = (0, 10),
                   # terminate_steady_state=true,
                   trip_lines = :dynamic,
                   trip_nodes = :dynamic,
                   trip_load_nodes = :none,
                   f_min = -2.5,
                   f_max = 1.5,
                   solverargs = (;dtmax=0.01), verbose = true);
    print("number of node failures "); print(length(sol.failures_nodes.saveval)); print("\n")
end
