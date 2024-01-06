# @assert VERSION == v"1.8.4"
const ON_YOGA = occursin("Yoga", gethostname())
const ON_PIK_HPC = occursin("login", gethostname())
const ON_POOL = occursin("pool", gethostname())

@info "Initialize environment"
# PKG_DIR = ON_YOGA ? abspath(@__DIR__, "..", "..", "..") : "/home/brandner/DynamicCascades.jl"

println("ON_PIK_HPC=$ON_PIK_HPC")

if ON_YOGA
    PKG_DIR = abspath(@__DIR__, "..", "..", "..")
    server_string = "YOGA_"
elseif ON_PIK_HPC
    PKG_DIR = "/home/brandner/DynamicCascades.jl"
    server_string = "PIK_HPC_"
elseif ON_POOL
    PKG_DIR = "/users/stud/brandner/MA/repos/DynamicCascades.jl"
    server_string = "POOL_"
end

using Pkg
Pkg.activate(PKG_DIR)


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
