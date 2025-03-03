const ON_CLUSTER = occursin("login", gethostname())

using Pkg
PKG_DIR = ON_CLUSTER ? "/home/brandner/DynamicCascades.jl" : joinpath(@__DIR__, "..", "..", "..")
Pkg.activate(PKG_DIR)
# Pkg.instantiate()

using LinearAlgebra
print("Number of threads before setting"); print(LinearAlgebra.BLAS.get_num_threads()); print("\n")
BLAS.set_num_threads(1)
print("Number of threads after setting"); print(LinearAlgebra.BLAS.get_num_threads()); print("\n")

using Revise
using DynamicCascades
using NetworkDynamics
using Graphs
using MetaGraphs
using Unitful
using Statistics
using GraphMakie
using Colors
using DynamicCascades: PLOT_DIR
using CairoMakie
using Dates
using DataFrames
using CSV

freq_bound = round(1.0/(2*π), digits=2)
# create folder
t=now()
datetime = Dates.format(t, "yyyymmdd_HHMMSS.s") # https://riptutorial.com/julia-lang/example/20476/current-time
folder = string("/",datetime,"inertia_vs_line+node_failures_f_bound=$freq_bound only11")
directory = string(RESULTS_DIR,folder)
mkpath(directory)

damping = 0.1u"s"
network0 = import_system(:rtsgmlc; damping, scale_inertia = 0.2, tconst = 0.01u"s")
nd, = nd_model(network0)
ω_state_idxs = idx_containing(nd, "ω")
gen_node_idxs = map(s -> parse(Int, String(s)[4:end]), nd.syms[ω_state_idxs])

scale_inertia_values = [0.2, 0.5, 0.8, 1.1, 1.4, 1.7, 2.0, 3.0, 4.0, 5.1, 6.1, 7.1, 8.0, 9.0, 10.0, 15.0, 21.0] # varying parameter

df_all_failures = DataFrame()
df_all_failures_nodes = DataFrame()
@time for scale_inertia in scale_inertia_values
    network = import_system(:rtsgmlc; damping, scale_inertia, tconst = 0.01u"s")
    x_static = steadystate(network)
    println("Scaling of inertia $scale_inertia \n steady state \n $x_static \n ")
    number_failures = Float64[]
    number_failures_nodes = Float64[]

    sol = simulate(network;
                   initial_fail = Int[11],
                   init_pert = :line,
                   tspan = (0, 500),
                   trip_lines = :none,
                   trip_nodes = :dynamic,
                   trip_load_nodes = :none,
                   f_min = -freq_bound,
                   f_max = freq_bound,
                   solverargs = (;dtmax=0.01), verbose = true);
                   # solverargs = (;dtmax=0.01, dtmin=0.0001, force_dtmin=true, save_everystep=false), verbose = true);

    push!(number_failures, length(sol.failures.saveval)-1) # `-1` as we don't want to count the initial failure
    push!(number_failures_nodes, length(sol.failures_nodes.saveval))

    df_all_failures[!, string(scale_inertia)] = number_failures
    df_all_failures_nodes[!, string(scale_inertia)] = number_failures_nodes
end

# write failures for each node/line to .csv
CSV.write(string(directory,"/all_failures.csv"), df_all_failures)
CSV.write(string(directory,"/all_failures_nodes.csv"), df_all_failures_nodes)
