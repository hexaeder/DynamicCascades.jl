"""
helpers functions
"""

using Serialization

###
### helper functions used by WS- and RTS-job array framework
###

## Watts-Strogatz
export get_network_args, import_system_wrapper, get_network_args_stripped, string_network_args, string_metagraph_args

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

#= #NOTE see `get_network_args`: Rather do `get_network_args` without units and use another function like get_network_args_with units
# This way units are added and removed again. =#
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

## RTS-GMCL testcase
export RTS_get_network_args, RTS_import_system_wrapper, RTS_get_network_args_stripped, RTS_string_network_args, RTS_string_metagraph_args

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


###
### pre- and postprocessing functions for WS- and RTS-job array framework
###
using Dates
export preprocess_WS, postprocess_jarray_data

#= TODO 
 - This is much slower with new ND. Not thoroughly tested why this is the case. Set `res_tol=1e-6` for small performance increase.
 - use kwargs!!!
=# 
function preprocess_WS(complement_to_existing_exp, existing_exp_name, name, exp_name, long_name,
    save_network_data, solver_name, steadystate_choice, N_ensemble_size, k_vals, β_vals, N_nodes, 
    inertia_values, K_vals, γ_vals, relate_inertia_and_damping, γ_eq_sq_I, γ_eq_I, τ_vals, σ_vals, μ_vals,
    failure_modes, gen_model, init_pert, freq_bounds, α_vals, monitored_power_flow;
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
        :save_network_data => save_network_data,
        :exp_name => exp_name, :long_name => long_name,
        :solver_name => solver_name, :steadystate_choice => steadystate_choice,
        :N_ensemble_size => N_ensemble_size,
        :k => k_vals, :β => β_vals, :N_nodes => N_nodes, # TODO change to :β_vals => β_vals? and so on...
        :inertia_values => inertia_values, :K => K_vals,  :γ => γ_vals, :τ => τ_vals,
        :σ => σ_vals, :μ => μ_vals,
        :failure_modes => failure_modes,
        :gen_model => gen_model,
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
    df_hpe[!, :filepath_steady_state] .= "<filepath>"; df_hpe[!, :filepath_graph] .= "<filepath>"; df_hpe[!, :filepath_power_injections] .= "<filepath>"

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
            - if steady state within tolerance exists (this automatically checks if θ ∈ [-π,π]).
            - if flows in steady state exceed rating by starting test simulation (depends on α).
            =#
            trial_counter = 1
            steady_state_check_approved = false
            while steady_state_check_approved == false
                try
                    N,k,β,graph_seed_,μ,σ,distr_seed_,K,_,_,_,_,_,_,_,_,_ = get_network_args(df_hpe, task_id)
                    # SteadyState mathematically does not depend on α, M, γ, τ
                    #= #NOTE NOTE NOTE NOTE NOTE NOTE NOTE NOTE NOTE NOTE NOTE NOTE NOTE NOTE NOTE NOTE NOTE NOTE NOTE NOTE NOTE NOTE NOTE NOTE NOTE NOTE NOTE NOTE NOTE NOTE NOTE NOTE NOTE
                    Here `α` is fixed, i.e. SteadyStates are generated with all flows < α.
                    For simulations with smaller values of default α=0.7, we get different
                    steady states and the ensemble is biased.=#
                    α = minimum(α_vals)
                    network = import_system(:wattsstrogatz; N=N, k=k, β=β, graph_seed=graph_seed, μ=μ, σ=σ, distr_seed=distr_seed, K=K, α=α, M=1u"s^2",  γ=1u"s", τ=1u"s")
                    println("task_id = $task_id: Test steady state of network with K=$K,N=$N,k=$k,β=$β,μ=$μ,σ=$σ,graph_seed=$graph_seed_,distr_seed=$distr_seed_...")
                    if steadystate_choice == :rootfind
                        x_static = steadystate(network; zeroidx=1, res_tol=1e-7, relax_init_guess = false) 
                    elseif steadystate_choice == :relaxation
                        x_static = steadystate(network; zeroidx=1, res_tol=1e-7, relax_init_guess = true) 
                    end

                    #= #NOTE This is a very conservative test. Alternative (not significantly faster): call `nd_model_and_CB!`,
                    then check if initial loads exceed rating (see `simulate`) =#
                    simulate(network;
                        gen_model=gen_model,
                        x_static=x_static,
                        initial_fail = [1],
                        failtime = 0.1,
                        tspan = (0, 0.10001),
                        verbose = false,
                        warn = false);

                    # If check is approved
                    trial_counter = 1 # reset trial_counter
                    steady_state_check_approved = true

                    # Save steady state and graph, create directories for steady states.
                    if save_network_data == true

                        graph_params = "graph_seed=$graph_seed_,distr_seed=$distr_seed_,k=$k,β=$β,ensemble_element=$ensemble_element"

                        # steady sate
                        steady_state_folder_path = joinpath(graph_combinations_path, "steady_states_graphs")
                        ispath(steady_state_folder_path) || mkdir(steady_state_folder_path)
                        filepath_SS = joinpath(steady_state_folder_path, string(graph_params,".csv"))
                        # Assign filepath to df
                        df_hpe[task_id:end,:filepath_steady_state] .= filepath_SS
                        # save steady state
                        CSV.write(filepath_SS, Dict(:SteadyState => x_static))

                        # graph
                        graph_folder_path = joinpath(graph_combinations_path, "graphs")
                        ispath(graph_folder_path) || mkdir(graph_folder_path)
                        filepath_g = joinpath(graph_folder_path, string(graph_params,"_graph"))
                        # Assign filepath to df
                        df_hpe[task_id:end,:filepath_graph] .= filepath_g
                        # save graph
                        savegraph(filepath_g, network.graph)

                        # power injections
                        power_injections_folder_path = joinpath(graph_combinations_path, "power_injections")
                        ispath(power_injections_folder_path) || mkdir(power_injections_folder_path)
                        filepath_PI = joinpath(power_injections_folder_path, string(graph_params,".csv"))
                        # Assign filepath to df
                        df_hpe[task_id:end,:filepath_power_injections] .= filepath_PI
                        # save power injections
                        Pmech = get_prop(network, 1:nv(network), :Pmech)
                        Pload = get_prop(network, 1:nv(network), :Pload)
                        CSV.write(filepath_PI, Dict(:Pmech => Pmech, :Pload => Pload))
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

"""
Calculates mean and standard error.

"""
# CHECK normalized sum of lines and nodes again.
function postprocess_jarray_data(exp_name_date)
    exp_data_dir = joinpath(RESULTS_DIR, exp_name_date)

    # load config file, and parameters
    df_config = DataFrame(CSV.File(joinpath(exp_data_dir, "config.csv")))
    exp_params_dict = Serialization.deserialize(joinpath(exp_data_dir, "exp.params"))

    N_ensemble_size = exp_params_dict[:N_ensemble_size]
    num_parameter_combinations = Int(length(df_config[!,:ArrayTaskID])/N_ensemble_size)

    df_avg_error = deepcopy(df_config)

    # Keep only the first N_rows rows
    df_avg_error = df_avg_error[1:num_parameter_combinations, :]

    # add columns to df
    # normalized
    df_avg_error[!, :ensemble_avg_norm_avg_line_failures] .= NaN; df_avg_error[!, :ensemble_avg_norm_avg_node_failures] .= NaN;
    df_avg_error[!, :ensemble_SE_norm_avg_line_failures] .= NaN; df_avg_error[!, :ensemble_SE_norm_avg_node_failures] .= NaN;
    df_avg_error[!, :ensemble_SE_norm_avg_node_plus_line_failures] .= NaN

    # not normalized
    df_avg_error[!, :ensemble_avg_line_failures] .= NaN; df_avg_error[!, :ensemble_avg_node_failures] .= NaN;
    df_avg_error[!, :ensemble_SE_line_failures] .= NaN; df_avg_error[!, :ensemble_SE_node_failures] .= NaN;
    df_avg_error[!, :ensemble_SE_avg_node_plus_line_failures] .= NaN

    # DataFrame with failures for ALL ArrayTaskIDs
    df_all_failures = deepcopy(df_config)
    df_all_failures[!, :avg_line_failures] .= NaN ; df_all_failures[!, :avg_node_failures] .= NaN;
    df_all_failures[!, :norm_avg_line_failures] .= NaN ; df_all_failures[!, :norm_avg_node_failures] .= NaN;

    # Delete columns of `df_avg_error` and generate `network`
    #= #NOTE / ENHANCEMENT: [2024-02-04 So]
    Here, `network` is only generated once. This is possible as all networks in one
    experiment have the same number of lines and nodes. If other networks are used, where
    in an ensemble of networks the number of lines and nodes varies (e.g. `:erdosrenyi`),
    `network` has to be generated for each element of the ensemble. One approach
    would the be to save `ne(network)` and `nv(network)` in df_config while preprocessing
    (Straighforward to implement but not necessarily needed).
    =#
    if exp_name_date[1:2] == "WS"
        select!(df_avg_error, Not([:graph_seed, :distr_seed, :filepath_steady_state, :filepath_graph, :filepath_power_injections, :ensemble_element]))
        network = import_system_wrapper(df_config, 1)
    elseif exp_name_date[1:3] == "RTS"
        select!(df_avg_error, Not([:ensemble_element]))
        network = RTS_import_system_wrapper(df_config, 1)
    else
        error("Please choose existing `exp_name_date`.")
    end

    # Find numer of (potentially failing) generator nodes.
    nw = nd_model_and_CB!(network)
    idxs_init_swing = map(idx -> idx.compidx, vidxs(nw, :, "ω"))
    nr_gen_nodes = length(idxs_init_swing)

    for task_id in df_avg_error.ArrayTaskID
        # loop over all elements of an ensemble
        avg_line_failures_ensemble = Float64[]; avg_node_failures_ensemble = Float64[]
        norm_avg_line_failures_ensemble = Float64[]; norm_avg_node_failures_ensemble = Float64[]
        for i in 0:num_parameter_combinations:(length(df_config[!,:ArrayTaskID]) - 1)
            # try...catch is for execution of postprocessing while not all jobs have finished
            try
                exp_data = joinpath(RESULTS_DIR, exp_name_date)
                if exp_name_date[1:2] == "WS"
                    N,k,β,graph_seed,μ,σ,distr_seed,K,α,M,γ,τ,freq_bound,trip_lines,trip_nodes,init_pert,ensemble_element = get_network_args_stripped(df_config, (task_id + i))
                    # _,_,_,_,_,_,_,_,_,_,_,_,freq_bound,trip_lines,trip_nodes,_,_ = get_network_args_stripped(df_config, (task_id + i))
                    graph_combinations_path = joinpath(exp_data, "k=$k,β=$β")
                elseif exp_name_date[1:3] == "RTS"
                    _,_,_,freq_bound,trip_lines,trip_nodes,_,_ = RTS_get_network_args_stripped(df_config, (task_id + i))
                    graph_combinations_path = exp_data
                end


                failure_mode_string = joinpath(graph_combinations_path, "trip_lines=$trip_lines,trip_nodes=$trip_nodes")
                failure_mode_frequ_bound = joinpath(failure_mode_string, "trip_lines=$trip_lines,trip_nodes=$trip_nodes,freq_bound=$freq_bound")

                # NOTE USE THIS ONLY FOR SINGLE EXPERIMENT (i.e. not when merging two experiments) ###############
                if exp_name_date[1:2] == "WS"
                    filename = string("/", string_network_args(df_config, task_id + i), ".csv")
                elseif exp_name_date[1:3] == "RTS"
                    filename = string("/", RTS_string_network_args(df_config, task_id + i), ".csv")
                end

                df_result = DataFrame(CSV.File(string(failure_mode_frequ_bound, filename)))
                ################################################################

                number_failures_lines = df_result[!, :number_failures_lines]
                number_failures_nodes = df_result[!, :number_failures_nodes]

                # average of a single run (one set of parameters, averaged over the number of lines of the network)
                # normalized
                norm_avg_line_failures = mean(number_failures_lines)/(ne(network)-1)
                norm_avg_node_failures = mean(number_failures_nodes)/nr_gen_nodes
                push!(norm_avg_line_failures_ensemble, norm_avg_line_failures)
                push!(norm_avg_node_failures_ensemble, norm_avg_node_failures)
                # not normalized
                avg_line_failures = mean(number_failures_lines)
                avg_node_failures = mean(number_failures_nodes)
                push!(avg_line_failures_ensemble, avg_line_failures)
                push!(avg_node_failures_ensemble, avg_node_failures)

                # normalized
                df_all_failures[(task_id + i), :norm_avg_line_failures] = norm_avg_line_failures
                df_all_failures[(task_id + i), :norm_avg_node_failures] = norm_avg_node_failures
                # not normalized
                df_all_failures[(task_id + i), :avg_line_failures] = avg_line_failures
                df_all_failures[(task_id + i), :avg_node_failures] = avg_node_failures

            catch
                continue
            end
        end
        # Calculate ensemble_avg and ensemble_standard_error and write to df
        # normalized
        df_avg_error[task_id,:ensemble_avg_norm_avg_line_failures] = mean(norm_avg_line_failures_ensemble)
        df_avg_error[task_id,:ensemble_avg_norm_avg_node_failures] = mean(norm_avg_node_failures_ensemble)
        df_avg_error[task_id,:ensemble_SE_norm_avg_line_failures] = 1 / sqrt(N_ensemble_size) * std(norm_avg_line_failures_ensemble; corrected=true)
        df_avg_error[task_id,:ensemble_SE_norm_avg_node_failures] = 1 / sqrt(N_ensemble_size) * std(norm_avg_node_failures_ensemble; corrected=true)
        df_avg_error[task_id,:ensemble_SE_norm_avg_node_plus_line_failures] = 1 / sqrt(N_ensemble_size) * std((norm_avg_line_failures_ensemble + norm_avg_node_failures_ensemble); corrected=true)

        # not normalized
        df_avg_error[task_id,:ensemble_avg_line_failures] = mean(avg_line_failures_ensemble)
        df_avg_error[task_id,:ensemble_avg_node_failures] = mean(avg_node_failures_ensemble)
        df_avg_error[task_id,:ensemble_SE_line_failures] = 1 / sqrt(N_ensemble_size) * std(avg_line_failures_ensemble; corrected=true)
        df_avg_error[task_id,:ensemble_SE_node_failures] = 1 / sqrt(N_ensemble_size) * std(avg_node_failures_ensemble; corrected=true)
        df_avg_error[task_id,:ensemble_SE_avg_node_plus_line_failures] = 1 / sqrt(N_ensemble_size) * std((avg_line_failures_ensemble + avg_node_failures_ensemble); corrected=true)
    end

    CSV.write(joinpath(RESULTS_DIR, exp_name_date, "avg_error.csv"), df_avg_error)
    CSV.write(joinpath(RESULTS_DIR, exp_name_date, "all_failures.csv"), df_all_failures)

    return df_avg_error, df_all_failures
end


###
### MetaGraph utils
###

"""
    describe_nodes(g::MetaGraph; firstcols=Vector{String}())

Returns DataFrame with all the node meta data.

NOTE NetworkDynamics.describe_vertices exists! 
DynamicCascades.describe_nodes(network) returns the properties
of the MetaGraph. NetworkDynamics.describe_vertices(nw) returns the 
properties of the NetworkDynamics.jl object.
"""
function describe_nodes(g::MetaGraph; firstcols=Vector{String}())
    df = DataFrame(; n=Int[])
    for n in 1:nv(g)
        row = push!(props(g, n), :n=>n)
        push!(df, row, cols=:union)
    end
    firstcols = String.(firstcols) # convert to string if given as Symbol
    has_prop(g, :NodeProps) &&
        append!(firstcols, String.(get_prop(g, :NodeProps)))
    append!(firstcols, names(df)) |> unique!
    select!(df, firstcols)
end

"""
    describe_edges(g::MetaGraph; firstcols=Vector{String}())

Returns DataFrame with all the edge meta data.

NOTE NetworkDynamics.describe_edges exists!
DynamicCascades.describe_edges(network) returns the properties
of the MetaGraph. NetworkDynamics.describe_edges(nw) returns the
properties of the NetworkDynamics.jl object.
"""
function describe_edges(g::MetaGraph; firstcols=Vector{String}())
    df = DataFrame(; src=Int[], dst=Int[])
    for e in edges(g)
        row = push!(props(g, e), :src=>e.src, :dst=>e.dst)
        push!(df, row, cols=:union)
    end
    firstcols = String.(firstcols) # convert to string if given as Symbol
    has_prop(g, :EdgeProps) &&
        append!(firstcols, String.(get_prop(g, :EdgeProps)))
    append!(firstcols, names(df)) |> unique!
    select!(df, firstcols)
end


using Makie
using Makie.GeometryBasics

import MetaGraphs: set_prop!, get_prop, has_prop

KEY_ITER = Union{AbstractUnitRange,Vector,AbstractEdgeIter}
"""
    set_prop!(g, keys::Iterable, prop::Symbol, vals::Iterable)

Set same property `prop` with different values `vals` for differet identifiers `keys`.
"""
function set_prop!(g, keys::KEY_ITER, prop::Symbol, vals::Vector)
    length(keys) == length(vals) || throw(ArgumentError("keys and vals needs to be of same length!"))
    for (k, val) in zip(keys, vals)
        if !ismissing(val)
            set_prop!(g, k, prop, val)
        end
    end
end

set_prop!(g, keys::KEY_ITER, p::Symbol, val) = set_prop!(g, keys, p, [val for v in keys])

"""
    get_prop(g, keys::Iterable, prop::Symbol)

Get same property `prop` with different values `vals` for differet keys.
"""
function get_prop(g, keys::KEY_ITER, prop::Symbol)
    [has_prop(g, k, prop) ? get_prop(g, k, prop) : missing for k in keys]
end

has_prop(g, keys::KEY_ITER, prop::Symbol) = all(k -> has_prop(g, k, prop), keys)


###
### analyze sol-object
###

using DataFrames

export describe_failures, collect_failure_times

"""
    describe_failures(sol)

Extracts failure times for both nodes and lines from `sol` and returns a DataFrame with columns
  - `:network_element`  : "name" of the failed vertex or edge,
  - `:idx`              : integer index of that vertex/edge,
  - `:failure_time`     : time of failure.
Rows are sorted by `failure_time`. 

Note that in legacy branch `mwe_old_ND_maybe_plots` the failure times were 
saved using `SavingCallback`.
"""
function describe_failures(sol)
    nw = NetworkDynamics.extract_nw(sol)

    # collect failure times for nodes and for lines
    failure_times_swing_nodes = collect_failure_times(sol,
                                    map(idx -> idx.compidx, vidxs(nw, :, "ω")), # indices that are initially swing nodes
                                    i -> vidxs(nw, i, :node_swing_stat))

    failure_times_lines = collect_failure_times(sol, 
                                    1:ne(nw),
                                    i -> eidxs(nw, i, :line_stat))

    # pair failure times with `network_element` in df
    rows_nodes = build_rows(failure_times_swing_nodes, NetworkDynamics.describe_vertices(nw))
    rows_lines = build_rows(failure_times_lines, NetworkDynamics.describe_edges(nw))

    df = DataFrame(vcat(rows_nodes, rows_lines))
    sort!(df, :failure_time)

    return df
end

#= helper: returns sorted Vector of with line/node index and failure time.
`indices`: line/node indices to be checked whether a failure occured
# `idxs_fun(i)`: e.g. `i -> vidxs(nw, i, :node_swing_stat)` =#
function collect_failure_times(sol, indices, idxs_fun)
    failure_times = Tuple{Int, Float64}[]
    for i in indices
        # check if element i has failed by the end
        if sol(sol.t[end], idxs = idxs_fun(i))[1] != 1 # see `SwingDynLoad`
            # loop through all (saved) time steps t to find the first j where line_stat/node_swing_stat != 1
            for j in eachindex(sol.t)
                if sol(sol.t[j], idxs = idxs_fun(i))[1] != 1
                    push!(failure_times, (i, sol.t[j]))
                    break
                end
            end
        end
    end
    sort!(failure_times, by = x -> x[2]) # sort by second column of tuple
    return failure_times
end

#= helper: connects the failure time with the name of the network element (swing_dyn_load, line, etc.)
that is failing and produces a Vector of NamedTuples (network_element, idx, failure_time)
desc_df = NetworkDynamics.describe_vertices(nw) or NetworkDynamics.describe_edges(nw)
with columns :idx and :name. =#
function build_rows(failure_times::Vector{Tuple{Int,Float64}}, desc_df::DataFrame)
    rows = NamedTuple[]
    for (idx, t) in failure_times
        name = desc_df.name[desc_df.idx .== idx][1]
        push!(rows, (network_element = name, idx = idx, failure_time = t))
    end
    return rows
end


export get_slurm_id
"""
Helper for reading out the output of a simulation on HPC.
Jobs are grouped in job arrays with same inertia constant allowing to adapt the slurm 
flag `--qos`. Thus `task_id` does not correspond to the slurm job_id. This function
returns the slurm job_id in order to manually read out the output of a simulation with 
given `task_id` in `/output`

Example: `WS_k=4_exp12_vary_I_only_lines_and_nodes_change_Pmech_complement_f_b,idx=2-3318425_44-csn14c171.out`
Here job_array_index=2, slurm_id=44

"""
function get_slurm_id(exp_name_date, task_id)
    exp_params_dict = Serialization.deserialize(joinpath(RESULTS_DIR, exp_name_date, "exp.params"))
    inertia_values = exp_params_dict[:inertia_values]
    N_inertia = length(inertia_values)

    df_config = DataFrame(CSV.File(joinpath(RESULTS_DIR, exp_name_date, "config.csv")))
    job_array_index = findfirst(x -> x == df_config[task_id, :inertia_values], inertia_values)

    slurm_id = (task_id - job_array_index) / N_inertia + 1

    println("job_array_index=$job_array_index, slurm_id=$slurm_id")
    return job_array_index, slurm_id
end