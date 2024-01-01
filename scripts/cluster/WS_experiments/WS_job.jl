# TODO adapt and include `@assert VERSION == v"1.6.0-rc2"

const ON_YOGA = occursin("Yoga", gethostname())

@info "Initialize environment on main process"

PKG_DIR = ON_YOGA ? abspath(@__DIR__, "..", "..", "..") : "/home/brandner/DynamicCascades.jl"
using Pkg
Pkg.activate(PKG_DIR)

if ON_YOGA
    using Revise
else # if on PIK-HPC or Pool
    #= TODO execute this in seperate script that is executed before this script.
    Otherwise `Pkg.instantiate()` is executed for every job.=#
    Pkg.instantiate()
    # Pkg.precompile()
end

using LinearAlgebra
print("Number of threads before setting"); print(LinearAlgebra.BLAS.get_num_threads()); print("\n")
BLAS.set_num_threads(1)
print("Number of threads after setting"); print(LinearAlgebra.BLAS.get_num_threads()); print("\n")

using DynamicCascades
using NetworkDynamics
using Graphs
using MetaGraphs
using Unitful
using Statistics
using GraphMakie
using Colors
using DynamicCascades: PLOT_DIR # TODO Probably remove
using CairoMakie
using Dates
using DataFrames
using CSV


exp_name_date = "WS_testrun_N_G=2_20240101_183724.014"

# read in SLURM_ARRAY_TASK_ID from `ARGS`
task_id = parse(Int64, ARGS[1])
# task_id = 1

# load config file
df_config = DataFrame(CSV.File(joinpath(@__DIR__, exp_name_date, "config.csv")))

# read in parameters from df_config
inertia = df_config[task_id,:inertia_values]
freq_bound = df_config[task_id,:freq_bounds]
trip_lines = Symbol(eval(Meta.parse(df_config[task_id,:failure_modes]))[1])
trip_nodes = Symbol(eval(Meta.parse(df_config[task_id,:failure_modes]))[2])
β = df_config[task_id,:beta]
k = df_config[task_id,:k]
N = df_config[task_id,:N_nodes]
K = df_config[task_id,:K]
γ = df_config[task_id,:gamma]
τ = df_config[task_id,:tau]
α = df_config[task_id,:alpha]
init_pert = Symbol(df_config[task_id,:init_pert])
σ = df_config[task_id,:sigma]
μ = df_config[task_id,:mu]
graph_seed = df_config[task_id,:graph_seed]
distr_seed = df_config[task_id,:distr_seed]

# Alternative of loading graphs that have been generated during postprocessing
# filepath_graph = df_config[task_id,:filepath]
# loadgraph(filepath_graph,MGFormat())

# SIMULATION ###################################################################

#= NOTE change inertia value of graph in case of loading graphs from .lg-files.
=> use correct inertia value via set_prop!(network, 1:nv(network), :_M, inertia_values[i] * 1u"s^2") =#


network = import_system(:wattsstrogatz; N=N, k=k, β=β, M=(inertia * 1u"s^2"), graph_seed=graph_seed,
    distr_seed=distr_seed, K=K, α=α, γ=(γ * 1u"s"), τ=(τ * 1u"s"), σ=σ)

# Find numer of (potentially failing) generator nodes.
if [get_prop(network,i,:type) for i in 1:nv(network)] == [:gen for i in 1:nv(network)]
    # This is the case for WS networks where initially all nodes are swing equation nodes.
    nr_gen_nodes = nv(network)
else
    # This is the case for the RTS testcases where initially NOT all nodes are swing equation nodes.
    nd, = nd_model(network)
    ω_state_idxs = idx_containing(nd, "ω")
    gen_node_idxs = map(s -> parse(Int, String(s)[4:end]), nd.syms[ω_state_idxs])
    nr_gen_nodes = length(gen_node_idxs)
end


number_failures_lines = Float64[]
number_failures_nodes = Float64[]
x_static = steadystate(network; verbose=true)
# for i in 1:ne(network)
for i in 1:4
    sol = simulate(network;
                   x_static=x_static,
                   initial_fail = Int[i],
                   init_pert = :line,
                   tspan = (0, 10),
                   trip_lines = trip_lines,
                   trip_nodes = trip_nodes,
                   trip_load_nodes = :none,
                   f_min = -freq_bound,
                   f_max = freq_bound,
                   solverargs = (;dtmax=0.01),
                   verbose = true);
    push!(number_failures_lines, length(sol.failures.saveval)-1) # `-1` as we don't want to count the initial failure
    push!(number_failures_nodes, length(sol.failures_nodes.saveval))
end
df_failures = DataFrame()
df_failures[!, :number_failures_lines] = number_failures_lines
df_failures[!, :number_failures_nodes] = number_failures_nodes

# calculate normalized average of failures
# line failures
df_failures[!, :norm_avg_line_failures] .= NaN
df_failures[1, :norm_avg_line_failures] = mean(number_failures_lines)/(ne(network)-1)
# node failures
df_failures[!, :norm_avg_node_failures] .= NaN
df_failures[1, :norm_avg_node_failures] = mean(number_failures_nodes)/nr_gen_nodes

# Create directory for results
exp_data = joinpath(RESULTS_DIR, exp_name_date)
ispath(exp_data) || mkdir(exp_data)

graph_combinations_path = joinpath(exp_data, "k=$k,beta=$β")
ispath(graph_combinations_path) || mkdir(graph_combinations_path)

failure_mode_string = joinpath(graph_combinations_path, "trip_lines=$trip_lines,trip_nodes=$trip_nodes")
ispath(failure_mode_string) || mkdir(failure_mode_string)
failure_mode_frequ_bound = joinpath(failure_mode_string, "trip_lines=$trip_lines,trip_nodes=$trip_nodes,freq_bound=$freq_bound")
ispath(failure_mode_frequ_bound) || mkdir(failure_mode_frequ_bound)

# write results to file
t=now()
datetime = Dates.format(t, "yyyymmdd_HHMMSS.s")

# TODO Hier string_network_args() verwenden
filename = 
"/trip_lines=$trip_lines,trip_nodes=$trip_nodes,freq_bound=$freq_bound,N=$N,k=$k,β=$β,graph_seed=$graph_seed,μ=$μ,σ=$σ,distr_seed=$distr_seed,K=$K,α=$α,M=$inertia,γ=$γ,τ=$τ,init_pert=$init_pert.csv"
CSV.write(string(failure_mode_frequ_bound, filename), df_failures)


# HIER WEITER #########################################################################################################
#=
 - /failure_mode=$failure_mode das überall einheitlich
 - paar simulationen als beispiel generieren
 - postprocessing und plotting coden
    - so coden, dass man postprocessing immer zwischendurch ausführen kann.
    - evtl. über task_id loopen und falls file nicht vorhanden `continue`
    - In Plot anzahl der elemente des Ensembles
=#


# # PLOTTING #####################################################################
#
# # load data
# df_inertia_vs_failures = DataFrame(CSV.File(string(directory,"/inertia_vs_failures.csv")))
# df_inertia_vs_failures_nodes = DataFrame(CSV.File(string(directory,"/inertia_vs_failures_nodes.csv")))
#
# # plot data
# fig = Figure(fontsize = 30)
# Axis(fig[1, 1],
#     # title = L"Decreasing $G_{av}$ for one sample grid",
#     # titlesize = 30,
#     xlabel = "inertia M [s^2]",
#     # xlabelsize = 30,
#     ylabel = "normalized average of failures",
#     # ylabelsize = 30
# )
#
# x = df_inertia_vs_failures.scale_inertia_values
# y_lines = df_inertia_vs_failures.rel_failures
# y_nodes = df_inertia_vs_failures_nodes.rel_failures
#
#
# scatter!(x, y_lines, label = "lines")
# scatter!(x, y_nodes, label = "nodes")
# axislegend(position = :rb)
#
# CairoMakie.save(string(directory,"/WS_inertia_vs_number_line+node_failures.pdf"),fig)
