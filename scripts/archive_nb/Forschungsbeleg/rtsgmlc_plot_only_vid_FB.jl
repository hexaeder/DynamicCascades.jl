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

f_bound = 1.0
sol = simulate(network;
               initial_fail = Int[27],
               # failtime=0.01,
               init_pert = :line,
               # P_perturb = 1.0,
               tspan = (0, 1.0),
               terminate_steady_state=true,
               trip_lines = :dynamic,
               trip_nodes = :dynamic,
               trip_load_nodes = :none,
               # f_min = -1.0/(2π),
               # f_max = 1.0/(2π),
               # f_min = -2.5,
               # f_max = 1.5,
               f_min = -f_bound,
               f_max = f_bound,
               solverargs = (;dtmax=0.01), verbose = true);

# plot solution
function plotnetwork(fig, sol, t)
    ax = Axis(fig)
    hidedecorations!(ax), hidespines!(ax)
    xlims!(ax, -10, 7)
    ylims!(ax, -5, 7)
    gpargs = gparguments(sol, t;
                         colortype=:abssteady,
                         Δω=Observable(f_bound),
                         offlinewidth=3,
                         offlinecolor=colorant"lightgray",
                         ecolorscaling = Observable(1.0),
                         node_size=15,
                         show_labels=false)
    p = graphplot!(ax, network; gpargs...)
    return ax, p
end

# tobs = Observable(sol.sol.t[end])
tobs = Observable(0.0)
fig = Figure(resolution=(1000,1000))

# plot network
# fig[1,1] = Label(fig, @lift("t = "*repr(round($tobs,digits=2))*" s"), tellwidth=false)
fig[0,:], p = plotnetwork(fig, sol, tobs)

fig[1,:] = Colorbar(fig, height=25, vertical=false,
                     colormap=Makie.ColorScheme([colorant"yellow", colorant"red"]), label="line load relative to rating")

# Does atm [2023-08-13 So] not work in combination with plotting node failures (offlinecolor gray)
# fig[2,:] = Colorbar(fig, get_node_plot(p), height=25, vertical=false, label="node frequency deviation [Hz]")

fig[2,:] = Colorbar(fig, height=25, limits = (-f_bound,f_bound), vertical=false,
                     colormap=ColorSchemes.diverging_bkr_55_10_c35_n256, label="node frequency deviation [Hz]")

# create video
# tobs = Observable(0.0)
T = 20 #20
tmax = sol.sol.t[end] #0.15 # tmax = 3.5
# tmax = 2
tmin = 0.0
fps = 30 # 30,20,100
trange = range(tmin, tmax, length=Int(T * fps))

record(fig, joinpath(MA_DIR,"rts_example_cascade_FB.mp4"), trange; framerate=30) do time
    tobs[] = time
end
