"""
Watts-Strogatz-Network-Ensemble: Using job array framework.
Heatmap: x: frequency bound f_b, y: Inertia I, z: log(sum of line + node failres + 1)
"""

include(abspath(@__DIR__, "..", "helpers_jarray.jl"))

if ON_YOGA
    using Revise
else # if on PIK-HPC or Pool
    Pkg.instantiate()
    # Pkg.precompile()
end

using GraphMakie
using Colors
using CairoMakie

# plotting parameters
create_posprocessing_data = false # set to `false` for fast plotting
sum_lines_nodes = true
normalize = false
heatmap_logscale = true
opacity = 0.3
fontsize = labelsize = 24
# markers
markersize = 15

exp_name_date = "WS_k=4_exp02_PIK_HPC_K_=3,N_G=32_20240208_000237.814"
exp_data_dir = joinpath(RESULTS_DIR, exp_name_date)
left_out_frequencies = [0.005, 0.01,
    0.16, 0.17, 0.18, 0.19, 0.2, 0.21, 0.22, 0.23, 0.24, 0.25, 0.26, 0.27, 0.28, 0.29, 0.3, 0.8]
left_out_inertia_values = []
left_out_inertia_values = [20.0, 30.0]
left_out_β_values = []

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
    select!(df_avg_error, Not([:graph_seed, :distr_seed, :filepath, :ensemble_element]))

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
    network = import_system_wrapper(df_config, 1)
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
                N,k,β,graph_seed,μ,σ,distr_seed,K,α,M,γ,τ,freq_bound,trip_lines,trip_nodes,init_pert,ensemble_element = get_network_args_stripped(df_config, (task_id + i))

                exp_data = joinpath(RESULTS_DIR, exp_name_date)
                graph_combinations_path = joinpath(exp_data, "k=$k,β=$β")

                failure_mode_string = joinpath(graph_combinations_path, "trip_lines=$trip_lines,trip_nodes=$trip_nodes")
                failure_mode_frequ_bound = joinpath(failure_mode_string, "trip_lines=$trip_lines,trip_nodes=$trip_nodes,freq_bound=$freq_bound")

                # NOTE use this only WHEN MERGING TWO EXPERIMENTS ###############
                # This is a bit hacky...
                file_path = ""
                test_counter = 0
                # read out all parameters , exept of seeds, read out ensemble_element find file in folder that matches both
                N,k,β,graph_seed,μ,σ,distr_seed,K,α,M,γ,τ,freq_bound,trip_lines,trip_nodes,init_pert,ensemble_element = get_network_args_stripped(df_config, (task_id + i))
                string1 = "trip_lines=$trip_lines,trip_nodes=$trip_nodes,freq_bound=$freq_bound,N=$N,k=$k,β=$β"
                string2 = "μ=$μ,σ=$σ"
                string3 = "K=$K,α=$α,M=$M,γ=$γ,τ=$τ,init_pert=$init_pert,ensemble_element=$ensemble_element."
                strings_to_match = [string1, string2, string3]
                files_in_dir = readdir(failure_mode_frequ_bound)
                for file_name in files_in_dir
                    if all(x -> x == true, [occursin(s, file_name) for s in strings_to_match])
                        file_path = joinpath(failure_mode_frequ_bound, file_name)
                        test_counter += 1
                        if test_counter > 1
                            error("Two files with the same parameters!")
                        end
                    end
                end
                df_result = DataFrame(CSV.File(file_path))
                ################################################################

                # USE THIS ONLY FOR SINGLE EXPERIMENT ###############################
                # filename = string("/", string_network_args(df_config, task_id + i), ".csv")
                # df_result = DataFrame(CSV.File(string(failure_mode_frequ_bound, filename)))
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
freq_bounds = exp_params_dict[:freq_bounds]

filtered_freq_bounds = filter!(x->x ∉ left_out_frequencies, deepcopy(freq_bounds))
filtered_inertia_values = filter!(x->x ∉ left_out_inertia_values, deepcopy(inertia_values))
filtered_β_values = filter!(x->x ∉ left_out_β_values, deepcopy(exp_params_dict[:β]))

# Create figures depending on the modes (loop).
failure_modes = exp_params_dict[:failure_modes]
β_vals = exp_params_dict[:β]

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
    N,k,β,graph_seed,μ,σ,distr_seed,K,α,M,γ,τ,freq_bound,trip_lines,trip_nodes,init_pert,ensemble_element = get_network_args_stripped(df_config, task_id)
    # frequency bounds
    if freq_bound ∈ left_out_frequencies
        continue
    end
    # inertia values
    if M ∈ left_out_inertia_values
        continue
    end
    # β values
    if β ∈ left_out_β_values
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
        N,k,β,graph_seed,μ,σ,distr_seed,K,α,M,γ,τ,freq_bound,trip_lines,trip_nodes,init_pert,ensemble_element = get_network_args_stripped(df_config, task_id)

        if (trip_lines == :dynamic &&  trip_nodes == :dynamic)
            if sum_lines_nodes == true
                # heatmap
                append!(all_failures_heatmap, (y_lines + y_nodes))
                # optimal inertia vs. f_b
                push!(min_failures, minimum(y_lines + y_nodes))
                push!(opt_inertia, filtered_inertia_values[argmin(y_lines + y_nodes)])
            end
        end
    end
end
N,k,β,graph_seed,μ,σ,distr_seed,K,α,M,γ,τ,freq_bound,trip_lines,trip_nodes,init_pert,ensemble_element = get_network_args_stripped(df_config, 1)

k_str = string(exp_params_dict[:k])
filtered_freq_bounds_str = string(filtered_freq_bounds)
K_str = string(exp_params_dict[:K])


if length(filtered_freq_bounds) > 1
    # create heatmap
    fig_hm = Figure(fontsize = (fontsize-3))
    ax_hm = Axis(fig_hm[1, 1], xticklabelrotation=π/2,
        title = "",
        xlabel = "Frequency bound f_b [Hz]", # frequency bound f_b
        xlabelsize = (fontsize + 5),
        ylabel = L"Inertia I [$s^2$]", # inertia value associated with minimum of failures
        ylabelsize = (fontsize + 5),
    )

    xs = filtered_freq_bounds
    ys = filtered_inertia_values
    data = transpose(reshape(all_failures_heatmap, length(filtered_inertia_values), length(filtered_freq_bounds)))

    # hm = heatmap!(ax_hm, xs, ys, data, colormap = Reverse(:blues))
    # hm = heatmap!(ax_hm, xs, ys, data, colormap = Reverse(color_map))
    # hm = heatmap!(ax_hm, xs, ys, heatmap_logscale ? log10.(data.+1) : data, colormap = Reverse(color_map))

    # Colormaps
    # color_map = ColorSchemes.plasma
    # color_map = ColorSchemes.cividis
    # color_map = :cividis
    # color_map = :blues
    # color_map = :grays

    hm = heatmap!(ax_hm, xs, ys, heatmap_logscale ? log10.(data.+1) : data, colormap = :grays)
    # fig_hm, ax_hm, hm = heatmap!(ax_hm, xs, ys, data, colormap = :blues)
    # fig_hm, ax, hm = heatmap(xs, ys, data)
    # https://docs.makie.org/stable/reference/blocks/colorbar/

    # create minimal failures (optimal inertia) vs. frequency bound f_b
    for i in 1:length(filtered_freq_bounds)
        if i == 1
            scatter!(ax_hm, filtered_freq_bounds[i], opt_inertia[i], color = Makie.wong_colors()[1], label = L"$I_{min}$: Inertia value associated with minimum of failures", markersize = markersize)
        end
        scatter!(ax_hm, filtered_freq_bounds[i], opt_inertia[i], color = Makie.wong_colors()[1], markersize = markersize)
    end
    axislegend(ax_hm, position = :rt, labelsize=(labelsize-8))
    Colorbar(fig_hm[:, end+1], hm, label = normalize ? L"normalized $N_{fail}$" : (heatmap_logscale ? L"$\log(N_{fail}+1)$" : L"$N_{fail}$"))

    ax_hm.xticks = filtered_freq_bounds
    ax_hm.xlabelpadding = 15
    ax_hm.yticks = [1.0, 3.0, 5.0, 7.5, 10.0, 20.0, 30.0]

    CairoMakie.save(joinpath(MA_DIR, "WS", "heatmap_log=$heatmap_logscale,sumlinesnodes=$sum_lines_nodes,K=$K_str,k=$k_str,β=$filtered_β_values,f_b_left_out=[_],M_left_out=$left_out_inertia_values.pdf"),fig_hm)
    CairoMakie.save(joinpath(MA_DIR, "WS", "heatmap_log=$heatmap_logscale,sumlinesnodes=$sum_lines_nodes,K=$K_str,k=$k_str,β=$filtered_β_values,f_b_left_out=[_],M_left_out=$left_out_inertia_values.png"),fig_hm)
end
fig_hm
