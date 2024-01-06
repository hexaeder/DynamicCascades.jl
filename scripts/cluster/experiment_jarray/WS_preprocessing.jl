include("helpers_jarray.jl")

# PARAMETERS ###################################################################
# Experiment name
save_graph_and_filepath = false
solver_name = "Rodas4P()" # NOTE adapt!
steadystate_choice = :rootfind # :relaxation
exp_name = "WS_testrun_paramsK=3_pool"
long_name = "all inertia, all β, one bound at 10" # for providing more details
# Graph params #############
N_nodes = 100
# k = [4, 10]
k = [4]
β = [0.1, 0.5, 0.9]
# β = [0.1, 0.5]

# MetaGraph params ###############
# inertia_values = [0.2, 0.5, 0.7, 0.9, 1.1, 1.4, 1.7, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 15.0, 20.0]
# inertia_values = [0.2, 1.0, 2.1, 3.0, 4.0, 5.0, 7.2, 11.0, 15.0, 21.0]
inertia_values = [0.01, 0.2, 1.0, 5.0, 10.0, 15.0]
# inertia_values = [0.2, 0.7, 5.0]
# inertia_values = [0.2, 0.7]

K = 9 # coupling K
γ = 1 # damping swing equation nodes γ
τ = 1 # time constant τ
σ = 1 # standard deviation σ
μ = 0 # mean μ

N_ensemble_size = 2 # 100

# Cacading params ##############

init_pert = [:line] # initial perturbation set constant to an initial line failure
α = 0.7 # tuning parameter α, :rating = α*K
monitored_power_flow = :apparent

#= frequency bounds [narrow bounds, wide bounds] bounds.
The value in numerator of round(0.1/(2*π) is the angular frequency =#

freq_bounds = [round(i/(2*π), digits=2) for i in [0.01, 0.1, 0.5, 1.0, 5.0]]

# freq_bounds = [round(0.1/(2*π), digits=2), round(0.5/(2*π), digits=2), round(10.0/(2*π), digits=2)]
# freq_bounds = [round(0.1/(2*π), digits=2)]

# failure_modes = [trip_lines, trip_nodes]
failure_modes = [[:dynamic, :dynamic], [:dynamic, :none], [:none, :dynamic]]
# failure_modes = [[:dynamic, :dynamic], [:dynamic, :none]]
# failure_modes = [[:dynamic, :dynamic]]


################################################################################



# Create result directory
t=now()
datetime = Dates.format(t, "_yyyymmdd_HHMMSS.s")
exp_name_date = string(exp_name, "_N_G=$N_ensemble_size", datetime)
exp_data_dir = joinpath(RESULTS_DIR, exp_name_date)
ispath(exp_data_dir) || mkdir(exp_data_dir)

# Writing parameters to files
exp_params_dict = Dict(
    :save_graph_and_filepath => save_graph_and_filepath,
    :exp_name => exp_name, :long_name => long_name,
    :solver_name => solver_name, :steadystate_choice => steadystate_choice,
    :N_ensemble_size => N_ensemble_size,
    :k => k, :β => β, :N_nodes => N_nodes,
    :inertia_values => inertia_values, :K => K,  :γ => γ, :τ => τ,
    :σ => σ, :μ => μ,
    :failure_modes => failure_modes,
    :init_pert => init_pert, :freq_bounds => freq_bounds, :α => α, :monitored_power_flow => monitored_power_flow,
    )

CSV.write(joinpath(exp_data_dir, "exp_params.csv"), exp_params_dict, writeheader=true, header=["parameter", "value"])
Serialization.serialize(joinpath(exp_data_dir, "exp.params"), exp_params_dict)



################################################################################
#= Create hyperparameter: the order is chosen such, that with an increasing number
of finished jobs of the job array the size of the ensemble increases equally,
e.g. when half of all jobs are finished, for each ensemble half of the grids
shall be simulated. Parameters included are parameters that are potentially changed
in an (future) experiment.=#
hyperparam = collect(Iterators.product(inertia_values, freq_bounds, failure_modes, β, k,
    N_nodes, K, γ, τ, α, init_pert, σ, μ))[:]

# Repeat hyperparam N_ensemble_size times
hyperparam_ensemble = repeat(hyperparam, N_ensemble_size)
# println(hyperparameter[1])

# create dataframe (hpe stands for hyperparam_ensemble)
#= NOTE `:inertia_values` needs to be in the second column, otherwise the big while
loop below is corrupted =#
df_hpe = DataFrame(map(idx -> getindex.(hyperparam_ensemble, idx), eachindex(first(hyperparam_ensemble))),
    [:inertia_values, :freq_bounds, :failure_modes, :β, :k, :N_nodes, :K, :γ, :τ, :α, :init_pert, :σ, :μ])

# add "ArrayTaskID" as first column of df
df_hpe = hcat(DataFrame(ArrayTaskID = 1:length(hyperparam_ensemble)), df_hpe)

# add columns `graph_seed`, `distr_seed` and `filepath` to df
df_hpe[!, :graph_seed] .= 0; df_hpe[!, :distr_seed] .= 0; df_hpe[!, :filepath] .= "<filepath>"

# For each row/ArrayTaskID of df_hpe add element of ensemble.
df_hpe[!, :ensemble_element] = vcat([fill(i, length(hyperparam)) for i in 1:N_ensemble_size]...)


# GENERATION OF NETWORKS #######################################################
number_of_task_ids_between_graphs = length(inertia_values) * length(freq_bounds) * length(failure_modes)
graph_seed = 0; distr_seed = 0
# Loop over each ArrayTaskID:
for task_id in df_hpe.ArrayTaskID

    #= For every new configuration of β and k and for each new element of an
    ensemble, generate a new network. The order of parameters in `hyperparam` is
    relevant for this.=#
    if ((task_id-1) % number_of_task_ids_between_graphs) == 0
        global graph_seed += 1; global distr_seed += 1
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
                    # TODO Avoid if-branch in loop.
                    if steadystate_choice == :rootfind
                        steadystate(network) # "Old" way: leads to some errors, thus the `catch`-option below
                    elseif steadystate_choice == :relaxation
                        steadystate_relaxation(network) # "New" way, steady state more precise, less/no errors, probabyl slower
                    end
                end

                # If (in case of no error in previous for loop) a steady state exists for all inertia values
                trial_counter = 1 # reset trial_counter
                steady_state_for_all_inertia_values = true

            # If a there is at least one inertia value for which no steady state exists:
            catch
                trial_counter += 1
                if trial_counter > max_trials
                    error("Tried $max_trials different values for `graph_seed` and `distr_seed`. Exiting...")
                    #= ENHANCEMENT: Instead of an error and an interuption of the
                    program leave out this grid and print which grid was left out.=#
                end

                N,k_,β,graph_seed_,μ,σ,distr_seed_,K_,_,_,γ,τ,_,_,_,_,_ = get_network_args_stripped(df_hpe, task_id)
                M = ustrip(u"s^2", get_prop(network, 1, :_M))
                # Next line needs to be kept with the long sting, string_metagraph_args() can't be used as `M` changes.
                @warn "No static solution found: ArrayTaskID=$task_id, Parameters: N=$N,k=$k_,β=$β,graph_seed=$graph_seed_,μ=$μ,σ=$σ,distr_seed=$distr_seed_,K=$K_,M=$M,γ=$γ,τ=$τ."

                # generate new grid by increasing seeds by one
                global graph_seed += 1; global distr_seed += 1
                df_hpe[task_id,:graph_seed] = graph_seed; df_hpe[task_id,:distr_seed] = distr_seed
                string_args = string_metagraph_args(df_hpe, task_id)
                println("Generate new network with changed seeds for ArrayTaskID=$task_id with parameters $string_args")
                network = import_system_wrapper(df_hpe, task_id)
            end
        end
    end
    df_hpe[task_id,:graph_seed] = graph_seed; df_hpe[task_id,:distr_seed] = distr_seed

    N,k_,β,graph_seed_,μ,σ,distr_seed_,K_,α,M,γ,τ,freq_bound,trip_lines,trip_nodes,_,ensemble_element = get_network_args_stripped(df_hpe, task_id)

    # Create directories for results (preventing that different jobs try to create a directory at the same time)
    graph_combinations_path = joinpath(exp_data_dir, "k=$k_,β=$β")
    ispath(graph_combinations_path) || mkdir(graph_combinations_path)

    failure_mode_string = joinpath(graph_combinations_path, "trip_lines=$trip_lines,trip_nodes=$trip_nodes")
    ispath(failure_mode_string) || mkdir(failure_mode_string)
    failure_mode_frequ_bound = joinpath(failure_mode_string, "trip_lines=$trip_lines,trip_nodes=$trip_nodes,freq_bound=$freq_bound")
    ispath(failure_mode_frequ_bound) || mkdir(failure_mode_frequ_bound)

    # Create directories for graphs
    if save_graph_and_filepath == true

        graph_folder_path = joinpath(graph_combinations_path, "graphs")
        ispath(graph_folder_path) || mkdir(graph_folder_path)

        graph_params = "graph_seed=$graph_seed_,distr_seed=$distr_seed_,k=$k_,β=$β,ensemble_element=$ensemble_element"
        filepath = joinpath(graph_folder_path, string(graph_params,".lg"))

        # Assign filepath to df
        df_hpe[task_id,:filepath] = relpath(filepath)

        # save network
        savegraph(filepath, network)
    end
end

# Save to CSV
CSV.write(joinpath(exp_data_dir, "config.csv"), df_hpe)
CSV.write(joinpath(exp_data_dir, "exp_params.csv"))
