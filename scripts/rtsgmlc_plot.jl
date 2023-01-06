#=
# Exemplary dynamic cascade in RTS-GMLC grid
![animation](rtsgmlc.mp4)
=#
# load used packages

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

using GLMakie#jl
GLMakie.activate!()#jl

# plot solution
function plotnetwork(fig, sol, t)
    ax = Axis(fig)
    hidedecorations!(ax), hidespines!(ax)
    xlims!(ax, -10, 7)
    ylims!(ax, -5, 7)
    gpargs = gparguments(sol, t;
                         colortype=:abssteady,
                         Δω=Observable(3.0),
                         offlinewidth=5,
                         offlinecolor=colorant"lightgray",
                         ecolorscaling = Observable(1.0),
                         node_size=15)
    p = graphplot!(ax, network; gpargs...)
    return ax, p
end


damping = 0.1u"s"
network = import_system(:rtsgmlc; damping, tconst = 0.01u"s")

x_static=steadystate(network; verbose=true),

init = 27
sol = simulate(network;
               initial_fail = Int[init],
               tspan = (0, 100),
               terminate_steady_state=true,
               trip_lines = :dynamic,
               solverargs = (;dtmax=0.01), verbose = true);

sol.network
sol.initial_fail
sol.failtime
sol.trip_lines

sol.load_S
sol.load_P
sol.failures
sol.sol

tobs = Observable(0.0)
fig = Figure(resolution=(1300,700))

fig[1,1] = Label(fig, @lift("t = "*repr(round($tobs,digits=2))*" s"), tellwidth=false)
fig[2,1], p = plotnetwork(fig, sol, tobs)

fig[:,2] = ax = Axis(fig; xlabel="time t in s", ylabel="apparent power flow in p.u.", title="flow transients of failing lines")
for (i, l) in pairs(sol.failures.saveval)
    t, s = seriesforidx(sol.load_S, l)
    lines!(ax, t, s; label="flow on edge ($i)", linewidth=5)
end
for (i, l) in pairs(sol.failures.saveval)
    t, s = seriesforidx(sol.load_S, l)
    scatter!(ax, (t[end], s[end]); marker=:star5, markersize=35)
end
vlines!(ax, tobs; color=:black, linewidth=1)

T = 20
tmax = 0.75
fps = 30
trange = range(0.0, tmax, length=Int(T * fps))

record(fig, "rtsgmlc.mp4", trange; framerate=30) do time
    tobs[] = time
end
