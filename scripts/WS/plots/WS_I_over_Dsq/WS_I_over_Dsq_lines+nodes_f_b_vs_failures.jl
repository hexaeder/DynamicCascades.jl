"""
Watts-Strogatz-Network-Ensemble: Using job array framework. Transition that appears
when varying the frequency bounds. Line and node failures summed.
"""
include(abspath(@__DIR__, "..", "..", "..", "helpers_jarray.jl"))

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

exp_name_date = "WS_k=4_exp05_3_I_over_Dsq_lines_and_nodes_PIK_HPC_K_=3,N_G=32_20250124_123933.775"
exp_data_dir = joinpath(RESULTS_DIR, exp_name_date)
left_out_frequencies = [
    0.16, 0.17, 0.18, 0.19, 0.2, 0.21, 0.22, 0.23, 0.24, 0.25, 0.26, 0.27, 0.28, 0.29, 0.3, 0.5, 0.8]
left_out_β_values = []

################################################################################
###################### Calculate mean and standard error #######################
################################################################################
if create_posprocessing_data == true
    postprocess_jarray_data(exp_name_date)
end
df_config = DataFrame(CSV.File(joinpath(exp_data_dir, "config.csv")))
exp_params_dict = Serialization.deserialize(joinpath(exp_data_dir, "exp.params"))


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


    fig_lines_and_nodes = Figure(size=(800,600),fontsize = fontsize)
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
        scatterlines!(ax_lines_and_nodes, filtered_freq_bounds, (y_lines + y_nodes), markersize = markersize, color=:black, label = "I=$i,k=$k,β=$β")
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
    lines!(ax_lines_and_nodes, [NaN], [NaN]; label=L"Damping $D=\sqrt{I}$ [$s$]", color=:white)
    axislegend(ax_lines_and_nodes, position = :rt, labelsize=labelsize)

    N,k,β,graph_seed,μ,σ,distr_seed,K,α,M,γ,τ,freq_bound,trip_lines,trip_nodes,init_pert,ensemble_element = get_network_args_stripped(df_config, 1)
    k_str = string(exp_params_dict[:k])
    filtered_freq_bounds_str = string(filtered_freq_bounds)
    K_str = string(exp_params_dict[:K])

    CairoMakie.save(joinpath(MA_DIR, "WS", "fb_vs_failures", "WS_I_over_Dsq_f_b_vs_failures_lines+nodes_sumlinesnodes=$sum_lines_nodes,K=$K_str,k=$k_str,β=$filtered_β_values,I=$i.png"),fig_lines_and_nodes)
    CairoMakie.save(joinpath(MA_DIR, "WS", "fb_vs_failures", "WS_I_over_Dsq_f_b_vs_failures_lines+nodes_sumlinesnodes=$sum_lines_nodes,K=$K_str,k=$k_str,β=$filtered_β_values,I=$i.pdf"),fig_lines_and_nodes)
    # fig_lines_and_nodes
end
