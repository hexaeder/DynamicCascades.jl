include("helpers_jarray.jl")

if ON_YOGA
    using Revise
# else # if on PIK-HPC or Pool
#     Pkg.instantiate() # For job arrays this leads to https://discourse.julialang.org/t/stale-file-handle-error-when-submitting-job-array-on-slurm/70108
#     # Pkg.precompile()
end

# PARAMETERS ###################################################################
exp_name_date = ARGS[2]

# load config file, and parameters
df_config = DataFrame(CSV.File(joinpath(RESULTS_DIR, exp_name_date, "config.csv")))
exp_params_dict = Serialization.deserialize(joinpath(RESULTS_DIR, exp_name_date, "exp.params"))

# read in SLURM_ARRAY_TASK_ID from `ARGS`
task_id = parse(Int64, ARGS[1])

# using job_array indices for splitting experiment into multiple job arrays
if length(ARGS) == 3
    job_array_index = parse(Int64, ARGS[3])
    N_inertia = length(exp_params_dict[:inertia_values])
    task_id = job_array_index + (task_id -1) * N_inertia
end

if length(ARGS) == 4
    job_array_index = parse(Int64, ARGS[3])
    freq_bound_index = parse(Int64, ARGS[4])
    N_inertia = length(exp_params_dict[:inertia_values])
    N_freq_bounds = length(exp_params_dict[:freq_bounds])
    task_id = job_array_index + (task_id -1) * N_inertia * N_freq_bounds + N_inertia * (freq_bound_index - 1)
end


M,γ,τ,freq_bound,trip_lines,trip_nodes,init_pert,ensemble_element = RTS_get_network_args_stripped(df_config, task_id)
monitored_power_flow = exp_params_dict[:monitored_power_flow]
steadystate_choice = exp_params_dict[:steadystate_choice]

# SIMULATION ###################################################################
# read in parameters from df_config
network = RTS_import_system_wrapper(df_config, task_id)

number_failures_lines = Float64[]
number_failures_nodes = Float64[]
if steadystate_choice == :rootfind
    x_static = steadystate(network; verbose=true) # "Old" way: leads to some errors, thus the `catch`-option below
elseif steadystate_choice == :relaxation
    x_static = steadystate_relaxation(network; verbose=true) # "New" way, steady state more precise, less/no errors, probabyl slower
end
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

# Write results to file
exp_data_dir = joinpath(RESULTS_DIR, exp_name_date)
failure_mode_frequ_bound = joinpath(exp_data_dir, "trip_lines=$trip_lines,trip_nodes=$trip_nodes", "trip_lines=$trip_lines,trip_nodes=$trip_nodes,freq_bound=$freq_bound")
filename = string("/", RTS_string_network_args(df_config, task_id), ".csv")
CSV.write(string(failure_mode_frequ_bound, filename), df_failures)
