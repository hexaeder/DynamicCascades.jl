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
using DynamicCascades: PLOT_DIR # TODO Probably remove
using Dates
using DataFrames
using CSV

# PARAMETERS ###################################################################
# Experiment name
save_graph_and_filepath = false
exp_name = "WS_testrun_N_G=2_"
t=now()
datetime = Dates.format(t, "yyyymmdd_HHMMSS.s")
exp_path = joinpath(@__DIR__, string(exp_name, datetime))
ispath(exp_path) || mkdir(exp_path)

# k = [4, 10]
k = [4]
beta = [0.1, 0.5, 0.9]
# beta = [0.1, 0.5]

#= frequency bounds [narrow bounds, wide bounds] bounds.
The value in numerator of round(0.1/(2*π) is the angular frequency =#
freq_bounds = [round(0.1/(2*π), digits=2), round(0.5/(2*π), digits=2)]
# freq_bounds = [round(0.1/(2*π), digits=2)]

# failure_modes = [trip_lines, trip_nodes]
failure_modes = [[:dynamic, :dynamic], [:dynamic, :none], [:none, :dynamic]]
# failure_modes = [[:dynamic, :dynamic], [:dynamic, :none]]
# failure_modes = [[:dynamic, :dynamic]]

# inertia_values = [0.2, 0.5, 0.7, 0.9, 1.1, 1.4, 1.7, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 15.0, 20.0]
# inertia_values = [0.2, 1.0, 2.1, 3.0, 4.0, 5.0, 7.2, 11.0, 15.0, 21.0]
# inertia_values = [0.2, 1.0, 5.0, 10.0, 15.0, 20.0]
inertia_values = [0.2, 0.7]

# constant parameters
N_ensemble_size = 2 # 100
N_nodes = 100

K = 6 # coupling K
gamma = 1 # damping swing equation nodes γ
tau = 1 # time constant τ

alpha = 0.7 # tuning parameter α, :rating = α*K
init_pert = [:line] # initial perturbation set constant to an initial line failure

sigma = 1 # standard deviation σ
mu = 0 # mean μ

################################################################################
#= create hyperparameter: the order is chosen such, that with an increasing number
of finished jobs of the job array the size of the ensemble increases equally,
e.g. when half of all jobs are finished, for each ensemble half of the grids
shall be simulated=#
hyperparam = collect(Iterators.product(inertia_values, freq_bounds, failure_modes, beta, k,
    N_nodes, K, gamma, tau, alpha, init_pert, sigma, mu))[:]

# Repeat hyperparam N_ensemble_size times
hyperparam_ensemble = repeat(hyperparam, N_ensemble_size)
# println(hyperparameter[1])

# TODO check
# TODO maybe separate the two Symbols: [:dynamic, :dynamic]
# create dataframe (hpe stands for hyperparam_ensemble)
#= NOTE `:inertia_values` needs to be in the second column, otherwise the big while
loop below is corrupted =#
df_hpe = DataFrame(map(idx -> getindex.(hyperparam_ensemble, idx), eachindex(first(hyperparam_ensemble))),
    [:inertia_values, :freq_bounds, :failure_modes, :beta, :k, :N_nodes, :K, :gamma, :tau, :alpha, :init_pert, :sigma, :mu])

# add "ArrayTaskID" as first column of df
df_hpe = hcat(DataFrame(ArrayTaskID = 1:length(hyperparam_ensemble)), df_hpe)

# add columns `graph_seed`, `distr_seed` and `filepath` filled with `nan` to df
df_hpe[!, :graph_seed] .= 0; df_hpe[!, :distr_seed] .= 0; df_hpe[!, :filepath] .= "<filepath>"


function get_network_args(df_hpe::DataFrame, task_id::Int)
    N=df_hpe[task_id,:N_nodes]
    k=df_hpe[task_id,:k]
    β=df_hpe[task_id,:beta]
    graph_seed=df_hpe[task_id,:graph_seed]
    μ=df_hpe[task_id,:mu]
    σ=df_hpe[task_id,:sigma]
    distr_seed=df_hpe[task_id,:distr_seed]
    K=df_hpe[task_id,:K]
    α=df_hpe[task_id,:alpha]
    M=df_hpe[task_id,:inertia_values]*1u"s^2"
    γ=df_hpe[task_id,:gamma]*1u"s"
    τ=df_hpe[task_id,:tau]*1u"s"
    freq_bound=df_hpe[task_id,:freq_bounds]
    failure_mode=df_hpe[task_id,:failure_modes]
    init_pert=df_hpe[task_id,:init_pert]

    return N,k,β,graph_seed,μ,σ,distr_seed,K,α,M,γ,τ,freq_bound,failure_mode,init_pert
end

#= NOTE Could be written as `import_system_wrapper(df_hpe::DataFrame, task_id::Int)`
using multiple dispatch, however not useful here in a script.=#
# TODO Change function name?
function import_system_wrapper(df_hpe::DataFrame, task_id::Int)
    N,k,β,graph_seed,μ,σ,distr_seed,K,α,M,γ,τ,_,_,_ = get_network_args(df_hpe, task_id)
    return import_system(:wattsstrogatz; N=N, k=k, β=β, graph_seed=graph_seed,
        μ=μ, σ=σ, distr_seed=distr_seed, K=K, α=α, M=M, γ=γ, τ=τ)
end

function string_network_args(df_hpe::DataFrame, task_id::Int)
    N,k,β,graph_seed,μ,σ,distr_seed,K,α,M,γ,τ,freq_bound,failure_mode,init_pert = get_network_args(df_hpe, task_id)
    trip_lines = failure_mode[1]; trip_nodes = failure_mode[2]
    return "trip_lines=$trip_lines,trip_nodes=$trip_nodes,freq_bound=$freq_bound,N=$N,k=$k,β=$β,graph_seed=$graph_seed,μ=$μ,σ=$σ,distr_seed=$distr_seed,K=$K,α=$α,M=$M,γ=$γ,τ=$τ,init_pert=$init_pert"
end


function string_metagraph_args(df_hpe::DataFrame, task_id::Int)
    N,k,β,graph_seed,μ,σ,distr_seed,K,_,M,γ,τ,_,_,_ = get_network_args(df_hpe, task_id)
    return "N=$N,k=$k,β=$β,graph_seed=$graph_seed,μ=$μ,σ=$σ,distr_seed=$distr_seed,K=$K,M=$M,γ=$γ,τ=$τ"
end

# GENERATION OF NETWORKS
number_of_task_ids_between_graphs = length(inertia_values) * length(freq_bounds) * length(failure_modes)

graph_seed = 0; distr_seed = 0
# Loop over each ArrayTaskID:
for task_id in df_hpe.ArrayTaskID

    #= For every new configuration of β and k and for each new element of an
    ensemble, generate a new network. The order of parameters in `hyperparam` is
    relevant for this.=#
    if ((task_id-1) % number_of_task_ids_between_graphs) == 0
        graph_seed += 1; distr_seed += 1
        df_hpe[task_id,:graph_seed] = graph_seed; df_hpe[task_id,:distr_seed] = distr_seed
        string_args = string_metagraph_args(df_hpe, task_id)
        println("Generate new MetaGraph:ArrayTaskID=$task_id with parameters $string_args")
        network = import_system_wrapper(df_hpe, task_id)

        #= Check for every inertia value if steady state within tolerance exists.
        If for one inertia value no steady state exists, a new MetaGraph is generated.=#
        max_trials = 10000
        trial_counter = 1
        steady_state_for_all_inertia_values = false
        while steady_state_for_all_inertia_values == false
            try
                # Try all inertia values
                for i in inertia_values
                    set_prop!(network, 1:nv(network), :_M, i * 1u"s^2")
                    println("Test steady state with inertia I=$i...")
                    steadystate(network)
                end

                # If (in case of no error in previous for loop) a steady state exists for all inertia values
                trial_counter = 1 # reset trial_counter
                steady_state_for_all_inertia_values = true

            # If a there is at least one inertia value for which no steady state exists
            catch
                trial_counter += 1
                if trial_counter > max_trials
                    error("Tried $max_trials different values for `graph_seed` and `distr_seed`. Exiting...")
                end

                N,k,β,graph_seed,μ,σ,distr_seed,K,_,_,γ,τ,_,_,_ = get_network_args(df_hpe, task_id)
                M = ustrip(u"s^2", get_prop(network, 1, :_M))
                @warn "No static solution found: ArrayTaskID=$task_id with parameters N=$N,k=$k,β=$β,graph_seed=$graph_seed,μ=$μ,σ=$σ,distr_seed=$distr_seed,K=$K,M=$M,γ=$γ,τ=$τ."

                # generate new grid by increasing seeds by one
                graph_seed += 1; distr_seed += 1
                df_hpe[task_id,:graph_seed] = graph_seed; df_hpe[task_id,:distr_seed] = distr_seed
                string_args = string_network_args(df_hpe, task_id)
                println("Generate new network with changed seeds for ArrayTaskID=$task_id with parameters $string_args")
                network = import_system_wrapper(df_hpe, task_id)
            end
        end
    end
    df_hpe[task_id,:graph_seed] = graph_seed; df_hpe[task_id,:distr_seed] = distr_seed

    N,k,β,graph_seed,μ,σ,distr_seed,K,α,M,γ,τ,freq_bound,failure_mode,init_pert = get_network_args(df_hpe, task_id)

    # Create paths and directories
    if save_graph_and_filepath == true

        graph_combinations_path = joinpath(exp_path, "k=$k,beta=$β")
        ispath(graph_combinations_path) || mkdir(graph_combinations_path)

        graph_folder_path = joinpath(graph_combinations_path, "graphs")
        ispath(graph_folder_path) || mkdir(graph_folder_path)

        graph_params = "graph_seed=$graph_seed,distr_seed=$distr_seed,k=$k,beta=$β"
        filepath = joinpath(graph_folder_path, string(graph_params,".lg"))

        # Assign filepath to df
        df_hpe[task_id,:filepath] = relpath(filepath)

        # save network
        savegraph(filepath, network)
    end


    # Create directories for results
    exp_data = joinpath(RESULTS_DIR, string(exp_name, datetime))
    ispath(exp_data) || mkdir(exp_data)

    graph_combinations_path = joinpath(exp_data, "k=$k,beta=$β")
    ispath(graph_combinations_path) || mkdir(graph_combinations_path)

    trip_lines = df_hpe[task_id,:failure_modes][1]
    trip_nodes = df_hpe[task_id,:failure_modes][2]

    failure_mode_string = joinpath(graph_combinations_path, "trip_lines=$trip_lines,trip_nodes=$trip_nodes")
    ispath(failure_mode_string) || mkdir(failure_mode_string)
    failure_mode_frequ_bound = joinpath(failure_mode_string, "trip_lines=$trip_lines,trip_nodes=$trip_nodes,freq_bound=$freq_bound")
    ispath(failure_mode_frequ_bound) || mkdir(failure_mode_frequ_bound)
end

# Save to CSV
CSV.write(joinpath(@__DIR__, string(exp_path, "/config.csv")), df_hpe)
