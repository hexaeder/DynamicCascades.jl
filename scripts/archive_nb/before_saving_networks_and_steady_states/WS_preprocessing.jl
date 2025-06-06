include("helpers_jarray.jl")

if ON_YOGA
    using Revise
else # if on PIK-HPC or Pool
    Pkg.instantiate()
    # Pkg.precompile()
end

# PARAMETERS ###################################################################
########### Only for adding new simulations to existing experiment #############
complement_to_existing_exp = false
# existing experiment
existing_exp_name = "WS_k=4_exp02_PIK_HPC_K_=3,N_G=32_20240208_000237.814"
# number of jobs that already have been simulated
N_jobs_total_existing = 864
# number of frquency bounds that are added
N_new_freq_bounds = length([2, 3, 4, 5, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40])
################################################################################

# Experiment name
name = "WS_k=4_exp02_"
long_name = "Uebergang frequency bounds further values for f_b." # for providing more details
save_graph_and_filepath = false
solver_name = "Rodas4P()"
steadystate_choice = :rootfind # :relaxation

# Graph params #############
N_nodes = 100
# k = [4, 10]
k_vals = [4]
β_vals = [0.5]
# β_vals = [0.25]
# β_vals = [0.5]

# MetaGraph params ###############
# inertia_values = [0.2, 0.5, 0.7, 0.9, 1.1, 1.4, 1.7, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 15.0, 20.0]
# inertia_values = [0.2, 1.0, 2.1, 3.0, 4.0, 5.0, 7.2, 11.0, 15.0, 21.0]
# inertia_values = [0.2, 1.0, 5.0, 10.0, 15.0]
inertia_values = [0.2, 0.5, 1.0, 3.0, 5.0, 7.5, 10.0, 20.0, 30.0]
K_vals = 3 # coupling K
γ_vals = 1 # damping swing equation nodes γ
τ_vals = 1 # time constant τ
σ_vals = 1 # standard deviation σ
μ_vals = 0 # mean μ

N_ensemble_size = 32 # 100

# Cascading params ##############
init_pert = [:line] # initial perturbation set constant to an initial line failure
α_vals = 0.7 # tuning parameter α, :rating = α*K
monitored_power_flow = :active

#= frequency bounds [narrow bounds, wide bounds] bounds.
The value in numerator of round(0.1/(2*π) is the angular frequency =#
# freq_bounds = [round(i/(2*π), digits=4) for i in [0.01, 0.1, 0.5, 1.0, 5.0]]
# freq_bounds = [round(0.1/(2*π), digits=2), round(0.5/(2*π), digits=2), round(10.0/(2*π), digits=2)]
# freq_bounds = [round(0.1/(2*π), digits=2)]
# freq_bounds = [0.005, 0.0075, 0.01, 0.015]#, 0.02, 0.025, 0.03]
# freq_bounds = [0.005, 0.03]
# freq_bounds = [0.005, 0.02, 0.03, 0.04, 0.8]
# This is frequency not angular frequency
freq_bounds = [0.005, 0.010, 0.015, 0.020, 0.025, 0.030, 0.035, 0.040, 0.045, 0.050, 0.055, 0.060, 0.065, 0.070, 0.075, 0.080, 0.085, 0.090, 0.095, 0.100, 0.110, 0.120, 0.130, 0.140, 0.150, 0.160, 0.170, 0.180, 0.190, 0.200, 0.210, 0.220, 0.230, 0.240, 0.250, 0.260, 0.270, 0.280, 0.290, 0.300, 0.800]

# freq_bounds = [0.02, 0.04]

# failure_modes = [trip_lines, trip_nodes]
# failure_modes = [[:dynamic, :dynamic], [:dynamic, :none], [:none, :dynamic]]
# failure_modes = [[:dynamic, :dynamic], [:dynamic, :none]]
failure_modes = [[:dynamic, :dynamic]]

exp_name_params = "K_=$K_vals,N_G=$N_ensemble_size"
exp_name = string(name, server_string, exp_name_params)
# exp_name = string(name, "PIK_HPC_", exp_name_params)

################################################################################

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
hyperparam = collect(Iterators.product(inertia_values, freq_bounds, failure_modes, β_vals, k_vals,
    N_nodes, K_vals, γ_vals, τ_vals, α_vals, init_pert, σ_vals, μ_vals))[:]

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


N_jobs_total = nrow(df_hpe)
N_inertia = length(inertia_values)
# job_array_length = Int64(N_jobs_total/N_inertia)
job_array_length = complement_to_existing_exp ? Int64((N_jobs_total - N_jobs_total_existing)/(N_inertia * N_new_freq_bounds)) : Int64(N_jobs_total/N_inertia)


exp_name_date_dict = Dict(
    :name => name,
    :exp_name_date => exp_name_date,
    :job_array_length => job_array_length,
    :N_inertia => N_inertia,
    )

CSV.write("sbatch_dict_$name.csv", exp_name_date_dict, writeheader=false)


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

                N,k_,β,graph_seed_,μ,σ,distr_seed_,K,_,_,γ,τ,_,_,_,_,_ = get_network_args_stripped(df_hpe, task_id)
                M = ustrip(u"s^2", get_prop(network, 1, :_M))
                # Next line needs to be kept with the long sting, string_metagraph_args() can't be used as `M` changes.
                # NOTE NOTE NOTE NOTE NOTE NOTE NOTE NOTE NOTE NOTE NOTE NOTE NOTE NOTE NOTE NOTE NOTE NOTE NOTE NOTE NOTE NOTE task_id stimmt hier nicht, weil man oben über i looped.
                @warn "No static solution found: ArrayTaskID=$task_id, Parameters: N=$N,k=$k_,β=$β,graph_seed=$graph_seed_,μ=$μ,σ=$σ,distr_seed=$distr_seed_,K=$K,M=$M,γ=$γ,τ=$τ."

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

    # Create directories for results (preventing that different jobs try to create a directory at the same time)
    N,k,β,graph_seed_,μ,σ,distr_seed_,K,α,M,γ,τ,freq_bound,trip_lines,trip_nodes,_,ensemble_element = get_network_args_stripped(df_hpe, task_id)
    graph_combinations_path = joinpath(exp_data_dir, "k=$k,β=$β")
    ispath(graph_combinations_path) || mkdir(graph_combinations_path)

    failure_mode_string = joinpath(graph_combinations_path, "trip_lines=$trip_lines,trip_nodes=$trip_nodes")
    ispath(failure_mode_string) || mkdir(failure_mode_string)
    failure_mode_frequ_bound = joinpath(failure_mode_string, "trip_lines=$trip_lines,trip_nodes=$trip_nodes,freq_bound=$freq_bound")
    ispath(failure_mode_frequ_bound) || mkdir(failure_mode_frequ_bound)

    # Create directories for graphs
    if save_graph_and_filepath == true

        graph_folder_path = joinpath(graph_combinations_path, "graphs")
        ispath(graph_folder_path) || mkdir(graph_folder_path)

        graph_params = "graph_seed=$graph_seed_,distr_seed=$distr_seed_,k=$k,β=$β,ensemble_element=$ensemble_element"
        filepath = joinpath(graph_folder_path, string(graph_params,".lg"))

        # Assign filepath to df
        df_hpe[task_id,:filepath] = relpath(filepath)

        # save network
        savegraph(filepath, network)
    end
end

# Save to CSV
CSV.write(joinpath(exp_data_dir, "config.csv"), df_hpe)
