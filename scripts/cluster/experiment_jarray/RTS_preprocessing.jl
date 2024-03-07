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
name = "RTS_exp02_"
long_name = "Uebergang ... and even more values for f_b." # for providing more details
solver_name = "Rodas4P()"
steadystate_choice = :rootfind # :relaxation


# MetaGraph params ###############
# inertia_values = [0.2, 0.5, 1.0, 3.0, 5.0, 7.5, 10.0, 20.0, 30.0]
inertia_values = [0.2, 0.5, 0.8, 1.1, 1.4, 1.7, 2.0, 3.0, 4.0, 5.1, 6.1, 7.1, 8.0, 9.0, 10.0, 15.0, 21.0]
γ_vals = 0.1 # damping swing equation nodes γ
τ_vals = 0.01 # time constant τ
N_ensemble_size = 1 # 100

# Cascading params ##############
init_pert = [:line] # initial perturbation set constant to an initial line failure

monitored_power_flow = :reactive

# frequency bounds
# freq_bounds = [0.40, 0.42, 0.44, 0.46, 0.48, 0.49, 0.50, 0.51, 0.52, 0.53, 0.54, 0.55, 0.56, 0.57, 0.58, 0.59, 0.60, 0.61, 0.62, 0.63, 0.64, 0.66, 0.68, 0.70]
freq_bounds = [0.01, 0.08, 0.10, 0.12, 0.14, 0.16, 0.18, 0.20, 0.22, 0.24, 0.26, 0.28, 0.30, 0.32, 0.34, 0.36, 0.38, 0.72, 0.74, 0.8, 0.85, 0.90, 1.00, 1.2, 1.4, 1.6, 1.8, 2.0]
# failure_modes = [trip_lines, trip_nodes]
# failure_modes = [[:dynamic, :dynamic], [:dynamic, :none], [:none, :dynamic]]
# failure_modes = [[:dynamic, :dynamic], [:dynamic, :none]]
failure_modes = [[:dynamic, :dynamic]]

exp_name = string(name, server_string)

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
    :exp_name => exp_name, :long_name => long_name,
    :solver_name => solver_name, :steadystate_choice => steadystate_choice,
    :N_ensemble_size => N_ensemble_size,
    :inertia_values => inertia_values, :γ => γ_vals, :τ => τ_vals,
    :failure_modes => failure_modes,
    :init_pert => init_pert, :freq_bounds => freq_bounds, :monitored_power_flow => monitored_power_flow,
    )

CSV.write(joinpath(exp_data_dir, "exp_params.csv"), exp_params_dict, writeheader=true, header=["parameter", "value"])
Serialization.serialize(joinpath(exp_data_dir, "exp.params"), exp_params_dict)



################################################################################
#= Create hyperparameter: the order is chosen such, that with an increasing number
of finished jobs of the job array the size of the ensemble increases equally,
e.g. when half of all jobs are finished, for each ensemble half of the grids
shall be simulated. Parameters included are parameters that are potentially changed
in an (future) experiment.=#
hyperparam = collect(Iterators.product(inertia_values, freq_bounds, failure_modes,
    γ_vals, τ_vals, init_pert))[:]

# Repeat hyperparam N_ensemble_size times
hyperparam_ensemble = repeat(hyperparam, N_ensemble_size)
# println(hyperparameter[1])

# create dataframe (hpe stands for hyperparam_ensemble)
#= NOTE `:inertia_values` needs to be in the second column, otherwise the big while
loop below is corrupted =#
df_hpe = DataFrame(map(idx -> getindex.(hyperparam_ensemble, idx), eachindex(first(hyperparam_ensemble))),
    [:inertia_values, :freq_bounds, :failure_modes, :γ, :τ, :init_pert])

# add "ArrayTaskID" as first column of df
df_hpe = hcat(DataFrame(ArrayTaskID = 1:length(hyperparam_ensemble)), df_hpe)

# For each row/ArrayTaskID of df_hpe add element of ensemble.
df_hpe[!, :ensemble_element] = vcat([fill(i, length(hyperparam)) for i in 1:N_ensemble_size]...)

N_jobs_total = nrow(df_hpe)
N_inertia = length(inertia_values)
job_array_length = complement_to_existing_exp ? Int64((N_jobs_total - N_jobs_total_existing)/(N_inertia * N_new_freq_bounds)) : Int64(N_jobs_total/N_inertia)

exp_name_date_dict = Dict(
    :name => name,
    :exp_name_date => exp_name_date,
    :job_array_length => job_array_length,
    :N_inertia => N_inertia,
    )

CSV.write("sbatch_dict_$name.csv", exp_name_date_dict, writeheader=false)


# Generate directories
for task_id in df_hpe.ArrayTaskID
    M,γ,τ,freq_bound,trip_lines,trip_nodes,_,ensemble_element = RTS_get_network_args_stripped(df_hpe, task_id)

    # Create directories for results (preventing that different jobs try to create a directory at the same time)
    graph_combinations_path = exp_data_dir
    ispath(graph_combinations_path) || mkdir(graph_combinations_path)

    failure_mode_string = joinpath(graph_combinations_path, "trip_lines=$trip_lines,trip_nodes=$trip_nodes")
    ispath(failure_mode_string) || mkdir(failure_mode_string)
    failure_mode_frequ_bound = joinpath(failure_mode_string, "trip_lines=$trip_lines,trip_nodes=$trip_nodes,freq_bound=$freq_bound")
    ispath(failure_mode_frequ_bound) || mkdir(failure_mode_frequ_bound)

end

# Save to CSV
CSV.write(joinpath(exp_data_dir, "config.csv"), df_hpe)
