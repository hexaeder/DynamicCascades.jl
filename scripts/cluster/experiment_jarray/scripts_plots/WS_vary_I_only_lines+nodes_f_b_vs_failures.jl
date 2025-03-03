"""
Watts-Strogatz-Network-Ensemble: Using job array framework. Transition that appears
when varying the frequency bounds. Line and node failures summed.
"""
#  NOTE NOTE NOTE NOTE NOTE NOTE NOTE NOTE NOTE NOTE NOTE NOTE Check normalized sum of lines and nodes again.
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
create_posprocessing_data = true # set to `false` for fast plotting
sum_lines_nodes = false
normalize = true
opacity = 0.30
fontsize = labelsize = 26
# markers
marker = (:circle, ":circle")
markersize = 15

exp_name_date = "WS_k=4_exp04_vary_I_only_lines_and_nodes_PIK_HPC_K_=3,N_G=32_20250126_012357.344"
exp_data_dir = joinpath(RESULTS_DIR, exp_name_date)
left_out_frequencies = [
    0.16, 0.17, 0.18, 0.19, 0.2, 0.21, 0.22, 0.23, 0.24, 0.25, 0.26, 0.27, 0.28, 0.29, 0.3, 0.8]
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
    select!(df_avg_error, Not([:graph_seed, :distr_seed, :filepath_steady_state, :ensemble_element]))

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
################################## Plotting ####################################
################################################################################
inertia_values = exp_params_dict[:inertia_values]
freq_bounds = exp_params_dict[:freq_bounds]

for i in inertia_values
    filtered_freq_bounds = filter!(x->x ∉ left_out_frequencies, deepcopy(freq_bounds))
    filtered_β_values = filter!(x->x ∉ left_out_β_values, deepcopy(exp_params_dict[:β]))
    # filtered_inertia_values = filter!(x->x ∉ left_out_inertia_values, deepcopy(inertia_values))
    left_out_inertia_values = filter(x -> x != i, deepcopy(inertia_values))


    fig_lines_and_nodes = Figure(fontsize = fontsize)
    ax_lines_and_nodes = Axis(fig_lines_and_nodes[1, 1], xticklabelrotation=π/2,
    # title = sum_lines_nodes ? "Summed line and node failures" : "Line and node failures",
    xlabel = "Frequency bound f_b [Hz]",
    ylabel = normalize ? "Normalized average of failures" : L"Averaged failures $N_{fail}$",
    )

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
    for task_id in df_avg_error.ArrayTaskID
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
    end

    filtered_freq_bounds
    (y_lines + y_nodes)
    if sum_lines_nodes == true
        scatterlines!(ax_lines_and_nodes, filtered_freq_bounds, (y_lines + y_nodes), markersize = markersize, color=:black, label = "I=$i,k=$k,β=$β,D=1 [s]")
        band!(ax_lines_and_nodes, filtered_freq_bounds, y_lines + y_nodes + err_nodes_plus_lines, y_lines + y_nodes - err_nodes_plus_lines, color=(:black, opacity), transparency=true)
    else
        lines!(ax_lines_and_nodes, [NaN], [NaN]; label="I=$i,k=$k,β=$β", color=:white, linewidth=3)
        scatterlines!(ax_lines_and_nodes, filtered_freq_bounds, y_lines, linestyle=:dash, markersize = markersize, color=:black, label = "Line failures (dashed)")
        band!(ax_lines_and_nodes, filtered_freq_bounds, y_lines + err_lines, y_lines - err_lines, color=(:black, opacity), transparency=true)

        scatterlines!(ax_lines_and_nodes, filtered_freq_bounds, y_nodes, markersize = markersize, color=:purple, label = "Node failures (solid)")
        band!(ax_lines_and_nodes, filtered_freq_bounds, y_nodes + err_nodes, y_nodes - err_nodes, color=(:pink, opacity+0.2), transparency=true)
        # lines!(ax_lines_and_nodes, [NaN], [NaN]; label="node failures: ___ (solid) \nline failures:   ----- (dashed)", color=:white, linewidth=3)
        # text!(ax_lines_and_nodes, 5.0, 0.02, text = "node failures: solid lines ___ \n line failures: dashed lines -----", align = (:center, :center), textsize=25)
    end


    ax_lines_and_nodes.xticks = filtered_freq_bounds
    ax_lines_and_nodes.xlabelpadding = 15
        lines!(ax_lines_and_nodes, [NaN], [NaN]; label=L"Damping $D=1$ [$s$]", color=:white)
    axislegend(ax_lines_and_nodes, position = :rt, labelsize=labelsize)

    N,k,β,graph_seed,μ,σ,distr_seed,K,α,M,γ,τ,freq_bound,trip_lines,trip_nodes,init_pert,ensemble_element = get_network_args_stripped(df_config, 1)
    k_str = string(exp_params_dict[:k])
    filtered_freq_bounds_str = string(filtered_freq_bounds)
    K_str = string(exp_params_dict[:K])

    CairoMakie.save(joinpath(MA_DIR, "WS", "fb_vs_failures", "WS_vary_I_only_f_b_vs_failures_lines+nodes_sumlinesnodes=$sum_lines_nodes,K=$K_str,k=$k_str,β=$filtered_β_values,I=$i.png"),fig_lines_and_nodes)
    CairoMakie.save(joinpath(MA_DIR, "WS", "fb_vs_failures", "WS_vary_I_only_f_b_vs_failures_lines+nodes_sumlinesnodes=$sum_lines_nodes,K=$K_str,k=$k_str,β=$filtered_β_values,I=$i.pdf"),fig_lines_and_nodes)
    # fig_lines_and_nodes
end
