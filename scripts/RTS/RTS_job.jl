include(abspath(@__DIR__, "..", "helpers_jarray.jl"))

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


# SIMULATION ###################################################################
_,_,_,freq_bound,trip_lines,trip_nodes,_,_ = RTS_get_network_args_stripped(df_config, task_id)
gen_model = exp_params_dict[:gen_model]
solver = exp_params_dict[:solver]

# read in network (includes parameters) from df_config
network = RTS_import_system_wrapper(df_config, task_id)

# read in steady state. The steady state is the same for all simulations
x_static = DataFrame(CSV.File(abspath(@__DIR__, "..", "..", "data/RTS-GMLC/steady_state_RTS-GMLC_20250604.csv"))).SteadyState

number_failures_lines = Float64[]
number_failures_nodes = Float64[]
for i in 1:ne(network)
    sol = simulate(network;
                    gen_model=gen_model,
                    x_static=x_static,
                    initial_fail = Int[i],
                    tspan = (0, 100000),
                    trip_lines = trip_lines,
                    trip_nodes = trip_nodes,
                    freq_bound = freq_bound,
                    solver = solver,
                    solverargs = (;reltol=1e-8, abstol=1e-6),
                    verbose = true);

    nw = NetworkDynamics.extract_nw(sol)
    idxs_init_swing = map(idx -> idx.compidx, vidxs(nw, :, "ω")) # indices that are initially swing nodes
    all_failing_nodes_idxs = [i for i in idxs_init_swing if sol(sol.t[end], idxs=vidxs(i, :node_swing_stat))[1] != 1]
    all_failing_lines_idxs = [i for i in 1:ne(nw) if sol(sol.t[end], idxs=eidxs(i, :line_stat))[1] == 0]
    push!(number_failures_nodes, length(all_failing_nodes_idxs))
    push!(number_failures_lines, length(all_failing_lines_idxs)-1) # `-1` as we don't want to count the initial failure
end
df_failures = DataFrame()
df_failures[!, :number_failures_lines] = number_failures_lines
df_failures[!, :number_failures_nodes] = number_failures_nodes

# Write results to file
exp_data_dir = joinpath(RESULTS_DIR, exp_name_date)
failure_mode_frequ_bound = joinpath(exp_data_dir, "trip_lines=$trip_lines,trip_nodes=$trip_nodes", "trip_lines=$trip_lines,trip_nodes=$trip_nodes,freq_bound=$freq_bound")
filename = string("/", RTS_string_network_args(df_config, task_id), ".csv")
CSV.write(string(failure_mode_frequ_bound, filename), df_failures)
