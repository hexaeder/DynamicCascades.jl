"""
Helper functions used by WS- and RTS-job array framework.
"""

# @assert VERSION == v"1.8.4"
const ON_YOGA = occursin("L7440", gethostname())
const ON_PIK_HPC = occursin("cs", gethostname())
const ON_POOL = occursin("pool", gethostname())

@info "Initialize environment"
# PKG_DIR = ON_YOGA ? abspath(@__DIR__, "..", "..", "..") : "/home/brandner/DynamicCascades.jl"

if ON_YOGA
    PKG_DIR = abspath(@__DIR__, "..", "..", "..")
    server_string = "L7440_"
elseif ON_PIK_HPC
    PKG_DIR = "/home/brandner/DynamicCascades.jl"
    server_string = "PIK_HPC_"
elseif ON_POOL
    PKG_DIR = "/users/stud/brandner/MA/repos/DynamicCascades.jl"
    server_string = "POOL_"
end

using Pkg
Pkg.activate(PKG_DIR)

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

# Watts-Strogatz
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

# RTS-GMCL testcase
function RTS_get_network_args(df::DataFrame, task_id::Int)
    M=df[task_id,:inertia_values]*1u"s^2"
    γ=df[task_id,:γ]*1u"s"
    τ=df[task_id,:τ]*1u"s"
    freq_bound=df[task_id,:freq_bounds]
    trip_lines=Symbol(eval(Meta.parse(string(df[task_id,:failure_modes])))[1])
    trip_nodes=Symbol(eval(Meta.parse(string(df[task_id,:failure_modes])))[2])
    init_pert=Symbol(df[task_id,:init_pert])
    ensemble_element=df[task_id,:ensemble_element]

    return M,γ,τ,freq_bound,trip_lines,trip_nodes,init_pert,ensemble_element
end

function RTS_import_system_wrapper(df::DataFrame, task_id::Int)
    M,γ,τ,freq_bound,trip_lines,trip_nodes,init_pert,ensemble_element = RTS_get_network_args(df, task_id)
    M = ustrip(u"s^2", M)
    return import_system(:rtsgmlc; damping=γ, scale_inertia = M, tconst = τ)
end

# Removing units.
function RTS_get_network_args_stripped(df::DataFrame, task_id::Int)
    M_,γ_,τ_,freq_bound,trip_lines,trip_nodes,init_pert,ensemble_element = RTS_get_network_args(df, task_id)
    M = ustrip(u"s^2", M_)
    τ = ustrip(u"s", τ_)
    γ = ustrip(u"s", γ_)
    return M,γ,τ,freq_bound,trip_lines,trip_nodes,init_pert,ensemble_element
end

function RTS_string_network_args(df::DataFrame, task_id::Int)
    M,γ,τ,freq_bound,trip_lines,trip_nodes,init_pert,ensemble_element = RTS_get_network_args_stripped(df, task_id)
    return "trip_lines=$trip_lines,trip_nodes=$trip_nodes,freq_bound=$freq_bound,M=$M,γ=$γ,τ=$τ,init_pert=$init_pert,ensemble_element=$ensemble_element"
end

function RTS_string_metagraph_args(df::DataFrame, task_id::Int)
    M,γ,τ,_,_,_,_ = RTS_get_network_args_stripped(df, task_id)
    return "M=$M,γ=$γ,τ=$τ"
end

function preprocess(complement_to_existing_exp, existing_exp_name, exp_name, long_name,
    save_graph_and_filepath, solver_name, steadystate_choice, N_ensemble_size, k_vals, β_vals, N_nodes, 
    inertia_values, K_vals, γ_vals, τ_vals, σ_vals, μ_vals,
    failure_modes, node_failure_model, init_pert, freq_bounds, α_vals, monitored_power_flow;
    max_trials = 1000)

    # Create result directory
    t=now()
    datetime = Dates.format(t, "_yyyymmdd_HHMMSS.s")
    # exp_name_date = string(exp_name, datetime)
    exp_name_date = complement_to_existing_exp ? existing_exp_name : string(exp_name, datetime)

    exp_data_dir = joinpath(RESULTS_DIR, exp_name_date)
    ispath(exp_data_dir) || mkdir(exp_data_dir)

    # Create folder for output + error of cluster runs
    output_path = joinpath(exp_data_dir, "output")
    ispath(output_path) || mkdir(output_path)

    # Writing parameters to files
    exp_params_dict = Dict(
        :save_graph_and_filepath => save_graph_and_filepath,
        :exp_name => exp_name, :long_name => long_name,
        :solver_name => solver_name, :steadystate_choice => steadystate_choice,
        :N_ensemble_size => N_ensemble_size,
        :k => k_vals, :β => β_vals, :N_nodes => N_nodes, # TODO change to :β_vals => β_vals? and so on...
        :inertia_values => inertia_values, :K => K_vals,  :γ => γ_vals, :τ => τ_vals,
        :σ => σ_vals, :μ => μ_vals,
        :failure_modes => failure_modes,
        :node_failure_model => node_failure_model,
        :init_pert => init_pert, :freq_bounds => freq_bounds, :α => α_vals, :monitored_power_flow => monitored_power_flow,
        )

    CSV.write(joinpath(exp_data_dir, "exp_params.csv"), exp_params_dict, writeheader=true, header=["parameter", "value"])
    Serialization.serialize(joinpath(exp_data_dir, "exp.params"), exp_params_dict)

    ################################################################################
    #= Create hyperparameter: the order is chosen such, that with an increasing number
    of finished jobs of the job array the size of the ensemble increases equally,
    e.g. when half of all jobs are finished, for each ensemble half of the grids
    shall be simulated. Parameters included are parameters that are potentially changed
    in an (future) experiment.=#
    hyperparam = collect(Iterators.product(inertia_values, γ_vals, τ_vals, K_vals, α_vals,
        freq_bounds, failure_modes, init_pert, β_vals, k_vals, N_nodes, σ_vals, μ_vals))[:]

    # Repeat hyperparam N_ensemble_size times
    hyperparam_ensemble = repeat(hyperparam, N_ensemble_size)
    # println(hyperparameter[1])

    # create dataframe (hpe stands for hyperparam_ensemble)
    df_hpe = DataFrame(map(idx -> getindex.(hyperparam_ensemble, idx), eachindex(first(hyperparam_ensemble))),
        [:inertia_values, :γ, :τ, :K, :α, :freq_bounds, :failure_modes, :init_pert, :β, :k, :N_nodes,  :σ, :μ])

    # add "ArrayTaskID" as first column of df
    df_hpe = hcat(DataFrame(ArrayTaskID = 1:length(hyperparam_ensemble)), df_hpe)

    # add columns `graph_seed`, `distr_seed` and `filepath` to df
    df_hpe[!, :graph_seed] .= 0; df_hpe[!, :distr_seed] .= 0
    df_hpe[!, :filepath_steady_state] .= "<filepath>"; df_hpe[!, :filepath_graph] .= "<filepath>"

    # For each row/ArrayTaskID of df_hpe add element of ensemble.
    df_hpe[!, :ensemble_element] = vcat([fill(i, length(hyperparam)) for i in 1:N_ensemble_size]...)

    #= inertia_variation - relate damping and inertia
    parameters need to be changed after `df_hpe` is populated except of seeds.=#
    if relate_inertia_and_damping
        # γ = τ = √I
        if γ_eq_sq_I
            for task_id in df_hpe.ArrayTaskID
                df_hpe.γ[task_id] = round(sqrt(df_hpe.inertia_values[task_id]), digits=4)
                # df_hpe.τ[task_id] = round(sqrt(df_hpe.inertia_values[task_id]), digits=4)
            end
        end
        if γ_eq_I
            # γ = τ = I
            for task_id in df_hpe.ArrayTaskID
                df_hpe.γ[task_id] = round(df_hpe.inertia_values[task_id], digits=4)
                # df_hpe.τ[task_id] = round(df_hpe.inertia_values[task_id], digits=4)
            end
        end
    end

    N_jobs_total = nrow(df_hpe)
    N_inertia = length(inertia_values)
    # job_array_length = Int64(N_jobs_total/N_inertia)
    # NOTE: Needs to be adapted if `complement_to_existing_exp = true´
    job_array_length = complement_to_existing_exp ? Int64((N_jobs_total - N_jobs_total_existing)/(N_inertia * N_new_freq_bounds)) : Int64(N_jobs_total/N_inertia)

    # Needs to be done AFTER `df_hpe` is populated
    exp_name_date_dict = Dict(
        :name => name,
        :exp_name_date => exp_name_date,
        :job_array_length => job_array_length,
        :N_inertia => N_inertia,
        )

    CSV.write("sbatch_dict_$name.csv", exp_name_date_dict, writeheader=false)
    
    # GENERATION OF NETWORKS #######################################################
    number_of_task_ids_between_graphs = length(inertia_values)*length(γ_vals)*length(τ_vals)*length(α_vals)*length(freq_bounds)*length(failure_modes)*length(init_pert)
    x_static = Float64[]
    graph_seed = 0; distr_seed = 0
    # Loop over each ArrayTaskID:
    for task_id in df_hpe.ArrayTaskID
        # Create directories for results (preventing that different jobs try to create a directory at the same time)
        N,k,β,graph_seed_,μ,σ,distr_seed_,K,α,M,γ,τ,freq_bound,trip_lines,trip_nodes,_,ensemble_element = get_network_args_stripped(df_hpe, task_id)
        graph_combinations_path = joinpath(exp_data_dir, "k=$k,β=$β")
        ispath(graph_combinations_path) || mkdir(graph_combinations_path)

        failure_mode_string = joinpath(graph_combinations_path, "trip_lines=$trip_lines,trip_nodes=$trip_nodes")
        ispath(failure_mode_string) || mkdir(failure_mode_string)
        failure_mode_frequ_bound = joinpath(failure_mode_string, "trip_lines=$trip_lines,trip_nodes=$trip_nodes,freq_bound=$freq_bound")
        ispath(failure_mode_frequ_bound) || mkdir(failure_mode_frequ_bound)

        #= For every new configuration of β and k and for each new element of an
        ensemble, generate a new network. The order of parameters in `hyperparam` is
        relevant for this.=#
        if ((task_id-1) % number_of_task_ids_between_graphs) == 0
            #= Find last idx with these parameters (N,k,β,graph_seed_,μ,σ,distr_seed_,K)
            with `idx < task_id`. This is in order to make sure that we use
            the same `graph_seed` and `distr_seed` for different experiments where other
            parameters (e.g. D,I) are used but the graphs and steady states need to be
            the same.=#
            #=The whole thing works only if this is true. Otherwise the condition
            `(task_id - idx) > 1` is never fulfilled.=#
            if length(N_nodes)*length(k_vals)*length(β_vals)*length(μ_vals)*length(σ_vals)*length(K_vals) > 1
                N,k,β,_,μ,σ,_,K,_,_,_,_,_,_,_,_,_ = get_network_args(df_hpe, task_id)
                idx = findlast(idx ->
                    (idx < task_id &&
                    N==df_hpe[idx,:N_nodes] &&
                    k==df_hpe[idx,:k] &&
                    β==df_hpe[idx,:β] &&
                    μ==df_hpe[idx,:μ] &&
                    σ==df_hpe[idx,:σ] &&
                    K==df_hpe[idx,:K]),
                    df_hpe.ArrayTaskID)
                #= If `typeof(idx) == Nothing` then the seeds are not changed. In this
                case this set of parameters is not in the previous `task_id`s  and the seeds
                start at 1. E.g. this is always the case for `task_id=1`.=#
                if typeof(idx) != Nothing && (task_id - idx) > 1
                    graph_seed = df_hpe[idx, :graph_seed]; distr_seed = df_hpe[idx, :distr_seed]
                else
                    graph_seed = 0; distr_seed = 0
                end
            end

            graph_seed += 1; distr_seed += 1
            df_hpe[task_id:end, :graph_seed] .= graph_seed; df_hpe[task_id:end, :distr_seed] .= distr_seed
            string_args = string_metagraph_args(df_hpe, task_id)
            println("Generate new MetaGraph:ArrayTaskID=$task_id with parameters $string_args")

            #= Check...
            - if steady state within tolerance exists (this autmatically checks if θ ∈ [-π,π].
            - if flows in steady state exceed rating by starting test simulation (depends on α.
            =#
            trial_counter = 1
            steady_state_check_approved = false
            while steady_state_check_approved == false
                try
                    N,k,β,graph_seed_,μ,σ,distr_seed_,K,_,_,_,_,_,_,_,_,_ = get_network_args(df_hpe, task_id)
                    # SteadyState mathematically does not depend on α, M, γ, τ
                    #= NOTE NOTE NOTE NOTE NOTE NOTE NOTE NOTE NOTE NOTE NOTE NOTE NOTE NOTE NOTE NOTE NOTE NOTE NOTE NOTE NOTE NOTE NOTE NOTE NOTE NOTE NOTE NOTE NOTE NOTE NOTE NOTE NOTE
                    Here `α` is fixed, i.e. SteadyStates are generated with all flows < α.
                    For simulations with smaller values of default α=0.7, we get different
                    steady states and the ensemble is biased.=#
                    α = minimum(α_vals)
                    network = import_system(:wattsstrogatz; N=N, k=k, β=β, graph_seed=graph_seed, μ=μ, σ=σ, distr_seed=distr_seed, K=K, α=α, M=1u"s^2",  γ=1u"s", τ=1u"s")
                    println("task_id = $task_id: Test steady state of network with K=$K,N=$N,k=$k,β=$β,μ=$μ,σ=$σ,graph_seed=$graph_seed_,distr_seed=$distr_seed_...")
                    if steadystate_choice == :rootfind
                        x_static = steadystate(network; zeroidx=1)
                    elseif steadystate_choice == :relaxation
                        x_static = steadystate_relaxation(network; zeroidx=1)
                    end

                    simulate(network;
                        x_static=x_static,
                        initial_fail = [1],
                        tspan = (0, 0.2),
                        monitored_power_flow = monitored_power_flow,
                        solverargs = (;dtmax=0.01),
                        verbose = false,
                        warn = false);

                    # If check is approved
                    trial_counter = 1 # reset trial_counter
                    steady_state_check_approved = true

                    # Save steady state and graph, create directories for steady states.
                    if save_graph_and_filepath == true

                        # save steady
                        steady_state_folder_path = joinpath(graph_combinations_path, "steady_states_graphs")
                        ispath(steady_state_folder_path) || mkdir(steady_state_folder_path)

                        graph_params = "graph_seed=$graph_seed_,distr_seed=$distr_seed_,k=$k,β=$β,ensemble_element=$ensemble_element"
                        filepath_SS = joinpath(steady_state_folder_path, string(graph_params,".csv"))

                        # Assign filepath to df
                        df_hpe[task_id:end,:filepath_steady_state] .= filepath_SS

                        # save steady state
                        steady_state_dict = Dict(:SteadyState => x_static)
                        CSV.write(filepath_SS, steady_state_dict)

                        # save graph
                        graph_folder_path = joinpath(graph_combinations_path, "graphs")
                        ispath(graph_folder_path) || mkdir(graph_folder_path)

                        graph_params = "graph_seed=$graph_seed_,distr_seed=$distr_seed_,k=$k,β=$β,ensemble_element=$ensemble_element"
                        filepath_g = joinpath(graph_folder_path, string(graph_params,"_graph"))

                        # Assign filepath to df
                        df_hpe[task_id:end,:filepath_graph] .= filepath_g

                        # save graph
                        savegraph(filepath_g, network.graph)
                    end

                # If a there is at least one parameter configuration for which no steady state exists:
                catch
                    trial_counter += 1
                    if trial_counter > max_trials
                        error("Tried $max_trials different values for `graph_seed` and `distr_seed`. Exiting...")
                        #= ENHANCEMENT: Instead of an error and an interuption of the
                        program leave out this grid and print which grid was left out.=#
                    end
                    @warn "No static solution found!"

                    # generate new grid by increasing seeds by one
                    graph_seed += 1; distr_seed += 1
                    # df_hpe[task_id,:graph_seed] = graph_seed; df_hpe[task_id,:distr_seed] = distr_seed
                    df_hpe[task_id:end, :graph_seed] .= graph_seed; df_hpe[task_id:end, :distr_seed] .= distr_seed
                    println("Generate new network with seeds +1.")
                end
            end
        end
    end

    CSV.write(joinpath(exp_data_dir, "config.csv"), df_hpe)
end
