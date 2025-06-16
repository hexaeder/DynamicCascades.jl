include(abspath(@__DIR__, "..", "scripts/helpers_jarray.jl"))

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

# in in WS_master_experiment.sh indices der neuen f_b Werte übergeben
# evtl. auch übergeben, ob complementing run

N,k,β,graph_seed,μ,σ,distr_seed,K,α,M,γ,τ,freq_bound,trip_lines,trip_nodes,init_pert,ensemble_element = get_network_args_stripped(df_config, task_id)
monitored_power_flow = exp_params_dict[:monitored_power_flow]
steadystate_choice = exp_params_dict[:steadystate_choice]

# Alternative of loading graphs that have been generated during postprocessing
# filepath_graph = df_config[task_id,:filepath]
# network = loadgraph(filepath_graph, MGFormat())
#= NOTE Saving graphs takes really long during preprocessing and one needs to save
a graph for each parameter configuration.
NOTE hange inertia value of graph in case of loading graphs from .mg-files.
=> use correct inertia value via set_prop!(network, 1:nv(network), :_M, inertia_values[i] * 1u"s^2") =#

# SIMULATION ###################################################################
# read in network (includes parameters) from df_config
network = import_system_wrapper(df_config, task_id)

# read in steady state
steady_state_dict  = CSV.File(df_config[task_id,:filepath_steady_state])
x_static = steady_state_dict[:SteadyState]


number_failures_lines = Float64[]
number_failures_nodes = Float64[]
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
                   solverargs = (;dtmax=0.01), # CHECK this severly prolongs simulations
                   verbose = true);
    push!(number_failures_lines, length(sol.failures.saveval)-1) # `-1` as we don't want to count the initial failure
    push!(number_failures_nodes, length(sol.failures_nodes.saveval))
    # CHECK TODO Replace by
    # eindices = [i for i in 1:ne(network) if sol(sol.t[end], idxs=eidxs(i, :line_stat))[1] == 0]
    # vindices = [i for i in 1:nv(network) if sol(sol.t[end], idxs=vidxs(i, :node_swing_stat))[1] == 0]

end
df_failures = DataFrame()
df_failures[!, :number_failures_lines] = number_failures_lines
df_failures[!, :number_failures_nodes] = number_failures_nodes

# Write results to file
exp_data_dir = joinpath(RESULTS_DIR, exp_name_date)
failure_mode_frequ_bound = joinpath(exp_data_dir, "k=$k,β=$β", "trip_lines=$trip_lines,trip_nodes=$trip_nodes", "trip_lines=$trip_lines,trip_nodes=$trip_nodes,freq_bound=$freq_bound")
filename = string("/", string_network_args(df_config, task_id), ".csv")
CSV.write(string(failure_mode_frequ_bound, filename), df_failures)
