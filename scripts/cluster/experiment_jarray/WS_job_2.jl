const ON_YOGA = occursin("Yoga", gethostname())

@info "Initialize environment"
PKG_DIR = ON_YOGA ? abspath(@__DIR__, "..", "..", "..") : "/home/brandner/DynamicCascades.jl"

using Pkg
Pkg.activate(PKG_DIR)

if ON_YOGA
    using Revise
else # if on PIK-HPC or Pool
    Pkg.instantiate() # leads to https://discourse.julialang.org/t/stale-file-handle-error-when-submitting-job-array-on-slurm/70108
#     # Pkg.precompile()
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
using Dates
using DataFrames
using CSV
using Serialization

function get_network_args(df::DataFrame, task_id::Int)
    N=df[task_id,:N_nodes]
    k=df[task_id,:k]
    β=df[task_id,:β]
    graph_seed=df[task_id,:graph_seed]
    μ=df[task_id,:μ]
    σ=df[task_id,:σ]
    distr_seed=df[task_id,:distr_seed]
    K=df[task_id,:K]
    α=df[task_id,:α]
    M=df[task_id,:inertia_values]*1u"s^2"
    γ=df[task_id,:γ]*1u"s"
    τ=df[task_id,:τ]*1u"s"
    freq_bound=df[task_id,:freq_bounds]
    trip_lines=Symbol(eval(Meta.parse(string(df[task_id,:failure_modes])))[1])
    trip_nodes=Symbol(eval(Meta.parse(string(df[task_id,:failure_modes])))[2])
    init_pert=Symbol(df[task_id,:init_pert])
    ensemble_element=df[task_id,:ensemble_element]

    return N,k,β,graph_seed,μ,σ,distr_seed,K,α,M,γ,τ,freq_bound,trip_lines,trip_nodes,init_pert,ensemble_element
end

function import_system_wrapper(df::DataFrame, task_id::Int)
    N,k,β,graph_seed,μ,σ,distr_seed,K,α,M,γ,τ,_,_,_,_,_ = get_network_args(df, task_id)
    return import_system(:wattsstrogatz; N=N, k=k, β=β, graph_seed=graph_seed,
        μ=μ, σ=σ, distr_seed=distr_seed, K=K, α=α, M=M, γ=γ, τ=τ)
end

# Removing units.
function get_network_args_stripped(df::DataFrame, task_id::Int)
    N,k,β,graph_seed,μ,σ,distr_seed,K,α,M_,γ_,τ_,freq_bound,trip_lines,trip_nodes,init_pert,ensemble_element = get_network_args(df, task_id)
    M = ustrip(u"s^2", M_)
    τ = ustrip(u"s", τ_)
    γ = ustrip(u"s", γ_)
    return N,k,β,graph_seed,μ,σ,distr_seed,K,α,M,γ,τ,freq_bound,trip_lines,trip_nodes,init_pert,ensemble_element
end

function string_network_args(df::DataFrame, task_id::Int)
    N,k,β,graph_seed,μ,σ,distr_seed,K,α,M,γ,τ,freq_bound,trip_lines,trip_nodes,init_pert,ensemble_element = get_network_args_stripped(df, task_id)
    return "trip_lines=$trip_lines,trip_nodes=$trip_nodes,freq_bound=$freq_bound,N=$N,k=$k,β=$β,graph_seed=$graph_seed,μ=$μ,σ=$σ,distr_seed=$distr_seed,K=$K,α=$α,M=$M,γ=$γ,τ=$τ,init_pert=$init_pert,ensemble_element=$ensemble_element"
end

function string_metagraph_args(df::DataFrame, task_id::Int)
    N,k,β,graph_seed,μ,σ,distr_seed,K,_,M,γ,τ,_,_,_,_ = get_network_args_stripped(df, task_id)
    return "N=$N,k=$k,β=$β,graph_seed=$graph_seed,μ=$μ,σ=$σ,distr_seed=$distr_seed,K=$K,M=$M,γ=$γ,τ=$τ"
end


# PARAMETERS ###################################################################
exp_name_date = ARGS[2]

# read in SLURM_ARRAY_TASK_ID from `ARGS`
task_id = parse(Int64, ARGS[1])
# task_id = 1

# load config file, and parameters
df_config = DataFrame(CSV.File(joinpath(RESULTS_DIR, exp_name_date, "config.csv")))
exp_params_dict = Serialization.deserialize(joinpath(RESULTS_DIR, exp_name_date, "exp.params"))

N,k,β,graph_seed,μ,σ,distr_seed,K,α,M,γ,τ,freq_bound,trip_lines,trip_nodes,init_pert,ensemble_element = get_network_args_stripped(df_config, task_id)
monitored_power_flow = exp_params_dict[:monitored_power_flow]

# Alternative of loading graphs that have been generated during postprocessing
# filepath_graph = df_config[task_id,:filepath]
# loadgraph(filepath_graph,MGFormat())
#= NOTE change inertia value of graph in case of loading graphs from .lg-files.
=> use correct inertia value via set_prop!(network, 1:nv(network), :_M, inertia_values[i] * 1u"s^2") =#

# SIMULATION ###################################################################
# read in parameters from df_config
network = import_system_wrapper(df_config, task_id)

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
# x_static = steadystate(network; verbose=true)
x_static = steadystate_relaxation(network; verbose=true)
for i in 1:ne(network)
    sol = simulate(network;
                   x_static=x_static,
                   initial_fail = Int[i],
                   init_pert = init_pert,
                   tspan = (0, 100000),
                   trip_lines = trip_lines,
                   trip_nodes = trip_nodes,
                   trip_load_nodes = :none,
                   monitored_power_flow = monitored_power_flow,
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

# Write results to file
exp_data_dir = joinpath(RESULTS_DIR, exp_name_date)
failure_mode_frequ_bound = joinpath(exp_data_dir, "k=$k,β=$β", "trip_lines=$trip_lines,trip_nodes=$trip_nodes", "trip_lines=$trip_lines,trip_nodes=$trip_nodes,freq_bound=$freq_bound")
filename = string("/", string_network_args(df_config, task_id), ".csv")
CSV.write(string(failure_mode_frequ_bound, filename), df_failures)
