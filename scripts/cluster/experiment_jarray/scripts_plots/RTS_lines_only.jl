"""
RTS-GMLC-Testcase: Using job array framework. Plotting lines only.
"""

include(abspath(@__DIR__, "..", "helpers_jarray.jl"))

if ON_YOGA
    using Revise
else # if on PIK-HPC or Pool
    Pkg.instantiate()
    # Pkg.precompile()
end

using GraphMakie
using Colors, ColorSchemes
using CairoMakie

# plotting parameters

show_title = false
create_posprocessing_data = true # set to `false` for fast plotting
sum_lines_nodes = true
normalize = false
opacity = 0.3
fontsize = labelsize = 28
# markers
markersize = 15


exp_name_date = "RTS_exp03_PIK_HPC__20240307_224712.402"
exp_data_dir = joinpath(RESULTS_DIR, exp_name_date)
left_out_frequencies = []
left_out_inertia_values = []

################################################################################
###################### Calculate mean and standard error #######################
################################################################################
if create_posprocessing_data == true
    # load config file, and parameters
    df_config = DataFrame(CSV.File(joinpath(exp_data_dir, "config.csv")))
    exp_params_dict = Serialization.deserialize(joinpath(exp_data_dir, "exp.params"))

    N_ensemble_size = exp_params_dict[:N_ensemble_size]
    num_parameter_combinations = Int(length(df_config[!,:ArrayTaskID])/N_ensemble_size)

    df_avg_error = deepcopy(df_config)

    # Delete columns
    select!(df_avg_error, Not([:ensemble_element]))

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

    #= NOTE / ENHANCEMENT: [2024-02-04 So]
    Here, `network` is only generated once. This is possible as all WS networks in one
    experiment have the same number of lines and nodes. If other networks are used, where
    in an ensemble of networks the number of lines and nodes varies (e.g. `:erdosrenyi`),
    `network` has to be generated for each element of the ensemble. The optimal soulution
    would the be to save `ne(network)` and `nv(network)` in df_config while preprocessing
    (Straighforward to implement but not necessarily needed).
    =#
    network = RTS_import_system_wrapper(df_config, 1)
    # Find numer of (potentially failing) generator nodes.
    if [get_prop(network,i,:type) for i in 1:nv(network)] == [:gen for i in 1:nv(network)]
        # This is the case for WS networks where initially all nodes are swing equation nodes.
        nr_gen_nodes = nv(network)
    else
        # This is the case for the RTS testcases where initially NOT all nodes are swing equation nodes.
        nd, = nd_model(network)
        ω_state_idxs = idx_containing(nd, "ω")
        gen_node_idxs = map(s -> parse(Int, String(s)[4:end]), nd.syms[ω_state_idxs])
        nr_gen_nodes = length(gen_node_idxs)
    end


    for task_id in df_avg_error.ArrayTaskID
        # loop over all elements of an ensemble
        avg_line_failures_ensemble = Float64[]; avg_node_failures_ensemble = Float64[]
        norm_avg_line_failures_ensemble = Float64[]; norm_avg_node_failures_ensemble = Float64[]
        for i in 0:num_parameter_combinations:(length(df_config[!,:ArrayTaskID]) - 1)
            # try...catch is for execution of postprocessing while not all jobs have finished
            try
                M,γ,τ,freq_bound,trip_lines,trip_nodes,init_pert,ensemble_element = RTS_get_network_args_stripped(df_config, (task_id + i))

                exp_data = joinpath(RESULTS_DIR, exp_name_date)
                graph_combinations_path = exp_data

                failure_mode_string = joinpath(graph_combinations_path, "trip_lines=$trip_lines,trip_nodes=$trip_nodes")
                failure_mode_frequ_bound = joinpath(failure_mode_string, "trip_lines=$trip_lines,trip_nodes=$trip_nodes,freq_bound=$freq_bound")


                # USE THIS ONLY FOR SINGLE EXPERIMENT ###############################
                filename = string("/", RTS_string_network_args(df_config, task_id + i), ".csv")
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
end

################################################################################
################################ Plotting  #####################################
################################################################################
inertia_values = exp_params_dict[:inertia_values]
filtered_inertia_values = filter!(x->x ∉ left_out_inertia_values, deepcopy(inertia_values))

fig_lines_only = Figure(fontsize = fontsize)
ax_lines_only = Axis(fig_lines_only[1, 1],
    title = show_title ? "Line failures" : "",
    xlabel = "Scaling factor of inertia I",
    ylabel = L"Averaged line failures $N_{fail}^L$",
)


# Create figures depending on the modes (loop).
failure_modes = exp_params_dict[:failure_modes]

df_avg_error = DataFrame(CSV.File(joinpath(RESULTS_DIR, exp_name_date, "avg_error.csv")))
inertia_values = exp_params_dict[:inertia_values]


#= Different inertia values for: Ensemble average over normalized average of
failures (the latter for a single network) =#
y_lines = Float64[]; y_nodes = Float64[]
#= Different inertia values for: Ensemble standard error over normalized average
of failures (the latter for a single network) =#
err_lines = Float64[]; err_nodes = Float64[]; err_nodes_plus_lines = Float64[]
all_failures_heatmap = Float64[]
min_failures = Float64[]; opt_inertia = Float64[]
for task_id in df_avg_error.ArrayTaskID # TODO renane variables: this is not an ArrayTaskID in the strict sense but an average over task IDs
    #= Empty arrays (after all inertia values of one configuration pushed to array)
    The entries in df_avg_error are ordered accordningly.=#
    if ((task_id-1) % length(inertia_values)) == 0
        y_lines = Float64[]; y_nodes = Float64[]
        err_lines = Float64[]; err_nodes = Float64[]; err_nodes_plus_lines = Float64[]
    end

    # Leave certain jobs out:
    M,γ,τ,freq_bound,trip_lines,trip_nodes,init_pert,ensemble_element = RTS_get_network_args_stripped(df_config, task_id)
    # frequency bounds
    if freq_bound ∈ left_out_frequencies
        continue
    end
    # inertia values
    if M ∈ left_out_inertia_values
        continue
    end

    # Read out ensemble_avg and ensemble_standard_error
    if normalize == true
        push!(y_lines, df_avg_error[task_id, :ensemble_avg_norm_avg_line_failures])
        push!(y_nodes, df_avg_error[task_id, :ensemble_avg_norm_avg_node_failures])
        push!(err_lines, df_avg_error[task_id, :ensemble_SE_norm_avg_line_failures])
        push!(err_nodes, df_avg_error[task_id, :ensemble_SE_norm_avg_node_failures])
        push!(err_nodes_plus_lines, df_avg_error[task_id,:ensemble_SE_norm_avg_node_plus_line_failures])
    else
        push!(y_lines, df_avg_error[task_id, :ensemble_avg_line_failures])
        push!(y_nodes, df_avg_error[task_id, :ensemble_avg_node_failures])
        push!(err_lines, df_avg_error[task_id, :ensemble_SE_line_failures])
        push!(err_nodes, df_avg_error[task_id, :ensemble_SE_node_failures])
        push!(err_nodes_plus_lines, df_avg_error[task_id,:ensemble_SE_avg_node_plus_line_failures])
    end

    # Only plot if all inertia values are pushed to `y_lines`, `y_nodes`, `err_lines`, `err_nodes`
    if M == maximum(filtered_inertia_values)
        # frequency argument first for a nice order in the legend
        M,γ,τ,freq_bound,trip_lines,trip_nodes,init_pert,ensemble_element = RTS_get_network_args_stripped(df_config, task_id)

        if (trip_lines == :dynamic &&  trip_nodes == :none)
            if freq_bound == filtered_freq_bounds[1]
                scatterlines!(ax_lines_only, filtered_inertia_values, y_lines, label = "", color = Makie.wong_colors()[1], linewidth = 3.5)
            end
        end
    end
end
M,γ,τ,freq_bound,trip_lines,trip_nodes,init_pert,ensemble_element = RTS_get_network_args_stripped(df_config, 1)

CairoMakie.save(joinpath(MA_DIR, "RTS_lines_only_M_left_out=$left_out_inertia_values.png"),fig_lines_only)
CairoMakie.save(joinpath(MA_DIR, "RTS_lines_only_M_left_out=$left_out_inertia_values.pdf"),fig_lines_only)
fig_lines_only
