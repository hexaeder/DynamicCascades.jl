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
name = "WS_k=4_exp03_I_over_Dsq_test_"
long_name = "Damping_variation_I_over_Dsq_." # for providing more details
save_graph_and_filepath = true
solver_name = "Rodas4P()"
steadystate_choice = :rootfind # :relaxation

# Graph params #############
N_nodes = 100
k_vals = [4]
β_vals = [0.5]

# MetaGraph params ###############
K_vals = 3 # coupling K

# NOTE see below "inertia_variation - relate inertia and damping"
inertia_values = [0.2, 0.5]
relate_inertia_and_damping = true
γ_eq_sq_I = true
γ_eq_I = false
if relate_inertia_and_damping
    γ_vals = NaN # damping swing equation nodes γ
else
    γ_vals = 1 # damping swing equation nodes γ
end
τ_vals = 1 # time constant τ

# Distribution power injections
σ_vals = 1 # standard deviation σ
μ_vals = 0 # mean μ

N_ensemble_size = 2

# Cascading params ##############
init_pert = [:line] # initial perturbation set constant to an initial line failure
α_vals = 0.7 # tuning parameter α, :rating = α*K
monitored_power_flow = :apparent

#= frequency bounds [narrow bounds, wide bounds] bounds.
The value in numerator of round(0.1/(2*π) is the angular frequency =#
# freq_bounds = [round(i/(2*π), digits=4) for i in [0.01, 0.1, 0.5, 1.0, 5.0]]
# This is frequency not angular frequency, 0.110, 0.120, 0.130, 0.140, 0.150, 0.160, 0.170, 0.180, 0.190, 0.200, 0.210, 0.220, 0.230, 0.240, 0.250, 0.260, 0.270, 0.280, 0.290, 0.300, 0.800
freq_bounds = [0.050, 0.100]

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
hyperparam = collect(Iterators.product(inertia_values, γ_vals, τ_vals, K_vals, α_vals,
    freq_bounds, failure_modes, init_pert, β_vals, k_vals, N_nodes, σ_vals, μ_vals))[:]

# Repeat hyperparam N_ensemble_size times
hyperparam_ensemble = repeat(hyperparam, N_ensemble_size)
# println(hyperparameter[1])

# create dataframe (hpe stands for hyperparam_ensemble)
#= NOTE `:inertia_values` needs to be in the second column, otherwise the big while
loop below is corrupted =#
df_hpe = DataFrame(map(idx -> getindex.(hyperparam_ensemble, idx), eachindex(first(hyperparam_ensemble))),
    [:inertia_values, :γ, :τ, :K, :α, :freq_bounds, :failure_modes, :init_pert, :β, :k, :N_nodes,  :σ, :μ])

# add "ArrayTaskID" as first column of df
df_hpe = hcat(DataFrame(ArrayTaskID = 1:length(hyperparam_ensemble)), df_hpe)

# add columns `graph_seed`, `distr_seed` and `filepath` to df
df_hpe[!, :graph_seed] .= 0; df_hpe[!, :distr_seed] .= 0; df_hpe[!, :filepath_steady_state] .= "<filepath>"

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


exp_name_date_dict = Dict(
    :name => name,
    :exp_name_date => exp_name_date,
    :job_array_length => job_array_length,
    :N_inertia => N_inertia,
    )

CSV.write("sbatch_dict_$name.csv", exp_name_date_dict, writeheader=false)


# GENERATION OF NETWORKS #######################################################
number_of_task_ids_between_graphs = length(inertia_values)*length(γ_vals)*length(τ_vals)*length(K_vals)*length(α_vals)*length(freq_bounds)*length(failure_modes)*length(init_pert)
number_of_parameter_configurations_relevant_for_steady_state = length(inertia_values)*length(γ_vals)*length(τ_vals)*length(K_vals)
graph_seed = 0; distr_seed = 0
# Loop over each ArrayTaskID:
for task_id in df_hpe.ArrayTaskID
    #= For every new configuration of β and k and for each new element of an
    ensemble, generate a new network. The order of parameters in `hyperparam` is
    relevant for this.=#
    if ((task_id-1) % number_of_task_ids_between_graphs) == 0
        global graph_seed += 1; global distr_seed += 1
        # df_hpe[task_id,:graph_seed] = graph_seed; df_hpe[task_id,:distr_seed] = distr_seed
        df_hpe[task_id:end, :graph_seed] .= graph_seed; df_hpe[task_id:end, :distr_seed] .= distr_seed
        string_args = string_metagraph_args(df_hpe, task_id)
        println("Generate new MetaGraph:ArrayTaskID=$task_id with parameters $string_args")

        #= Check...
         - if steady state within tolerance exists (this autmatically checks if θ ∈ [-π,π]).
         - if flows in steady state exceed rating by starting test simulation.
        =#
        max_trials = 10000
        trial_counter = 1
        steady_state_check_approved = false
        while steady_state_check_approved == false
            try
                network = import_system_wrapper(df_hpe, task_id)
                _,_,_,graph_seed_,_,_,distr_seed_,K,_,M,γ,τ,_,_,_,_,ensemble_element = get_network_args_stripped(df_hpe, task_id)
                println("task_id = $task_id: Test steady state of network with inertia I=$M,γ=$γ,τ=$τ,K=$K (graph_seed=$graph_seed_, distr_seed=$distr_seed_)...")
                if steadystate_choice == :rootfind
                    x_static = steadystate(network; tol=1e-6, zeroidx=1)
                elseif steadystate_choice == :relaxation
                    x_static = steadystate_relaxation(network; zeroidx=1)
                end

                simulate(network;
                   x_static=x_static,
                   initial_fail = [1],
                   tspan = (0, 0.2),
                   monitored_power_flow = monitored_power_flow,
                   solverargs = (;dtmax=0.01),
                   verbose = true);

                # If check is approved
                trial_counter = 1 # reset trial_counter
                steady_state_check_approved = true

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
                global graph_seed += 1; global distr_seed += 1
                # df_hpe[task_id,:graph_seed] = graph_seed; df_hpe[task_id,:distr_seed] = distr_seed
                df_hpe[task_id:end, :graph_seed] .= graph_seed; df_hpe[task_id:end, :distr_seed] .= distr_seed
                println("Generate new network with seeds +1.")
            end
        end
    end

    # Create directories for results (preventing that different jobs try to create a directory at the same time)
    N,k,β,graph_seed_,μ,σ,distr_seed_,K,α,M,γ,τ,freq_bound,trip_lines,trip_nodes,_,ensemble_element = get_network_args_stripped(df_hpe, task_id)
    graph_combinations_path = joinpath(exp_data_dir, "k=$k,β=$β")
    ispath(graph_combinations_path) || mkdir(graph_combinations_path)

    failure_mode_string = joinpath(graph_combinations_path, "trip_lines=$trip_lines,trip_nodes=$trip_nodes")
    ispath(failure_mode_string) || mkdir(failure_mode_string)
    failure_mode_frequ_bound = joinpath(failure_mode_string, "trip_lines=$trip_lines,trip_nodes=$trip_nodes,freq_bound=$freq_bound")
    ispath(failure_mode_frequ_bound) || mkdir(failure_mode_frequ_bound)

    # Save steady state, create directories for graphs, save graphs.
    if save_graph_and_filepath == true
        # save steady state
        set_prop!(network, :SteadyState, x_static)

        steady_state_folder_path = joinpath(graph_combinations_path, "steady_states_graphs")
        ispath(steady_state_folder_path) || mkdir(steady_state_folder_path)

        graph_params = "graph_seed=$graph_seed_,distr_seed=$distr_seed_,k=$k,β=$β,ensemble_element=$ensemble_element"
        filepath = joinpath(steady_state_folder_path, string(graph_params,".csv"))

        # Assign filepath to df
        df_hpe[task_id:end,:filepath_steady_state] .= filepath

        # # save network (this is really slow!)
        # savegraph(filepath, network)

        # save steady state
        steady_state_dict = Dict(:SteadyState => x_static)
        CSV.write(filepath, steady_state_dict)

    end
end

# Save to CSV
CSV.write(joinpath(exp_data_dir, "config.csv"), df_hpe)
