"""
Model: line and node failures.
Inertia vs. line and node failures summed for 3 different frequency bounds (narrow and intermediate bound).
"""

using GraphMakie
using Colors, ColorSchemes
using CairoMakie
CairoMakie.activate!()

# plotting parameters
create_posprocessing_data = false # set to `false` for fast plotting
sum_lines_nodes = true
normalize = false
opacity = 0.3
fontsize = labelsize = 40
# line_colors = [Makie.wong_colors()[3], Makie.wong_colors()[5], Makie.wong_colors()[6]]
line_colors = ["#5C3D99FF", "#56B4E9FF", "#D55E00FF"] # [deep indigo, Hellblau, Orange]
# markers
markersize = 15
linewidth = 4

exp_name_date = "RTS_exp04_variation_frequency+inertia_PIK_HPC_20250616_213442.721"
exp_data_dir = joinpath(RESULTS_DIR, exp_name_date)
left_out_frequencies = [0.01, 0.08, 0.1, 0.12, 0.14, 0.16, 0.18, 0.22, 0.24, 0.26, 0.28,
    0.32, 0.34, 0.36, 0.38, 0.40, 0.42, 0.44, 0.46, 0.48, 0.49, 0.50, 0.51, 0.52, 0.53, 0.54, 0.55,
    0.56, 0.57, 0.58, 0.59, 0.60, 0.61, 0.62, 0.63, 0.64, 0.66, 0.68, 0.70, 0.72, 0.74,
    0.8, 0.85, 0.90, 1.2, 1.4, 1.6, 1.8, 2.0] # combined runs

# left_out_frequencies = [0.01, 0.1, 0.08, 0.12, 0.14, 0.16, 0.18, 0.20, 0.22, 0.24, 0.26,
#     0.32, 0.34, 0.36, 0.38, 0.40, 0.42, 0.44, 0.46, 0.48, 0.49, 0.50, 0.51, 0.52, 0.53, 0.54, 0.55,
#     0.56, 0.57, 0.58, 0.59, 0.60, 0.61, 0.62, 0.63, 0.64, 0.66, 0.68, 0.70, 0.72, 0.74,
#     0.8, 0.85, 0.90, 1.2, 1.4, 1.6, 1.8, 2.0] # combined runs

left_out_inertia_values = [0.2]

################################################################################
###################### Calculate mean and standard error #######################
################################################################################
if create_posprocessing_data == true
    postprocess_jarray_data(exp_name_date)
end
df_config = DataFrame(CSV.File(joinpath(exp_data_dir, "config.csv")))
exp_params_dict = Serialization.deserialize(joinpath(exp_data_dir, "exp.params"))

################################################################################
################################ Plotting  #####################################
################################################################################
inertia_values = exp_params_dict[:inertia_values]
freq_bounds = exp_params_dict[:freq_bounds]

filtered_freq_bounds = filter!(x->x ∉ left_out_frequencies, deepcopy(freq_bounds))
filtered_inertia_values = filter!(x->x ∉ left_out_inertia_values, deepcopy(inertia_values))

fig_lines_and_nodes = Figure(size=(800,600),fontsize = fontsize)
ax_lines_and_nodes = Axis(fig_lines_and_nodes[1, 1],
    title = "Node-line failure model (PG)",
    titlefont = :regular,
    xlabel = "Inertia scaling factor",
    ylabel = normalize ? "normalized average of failures" : L"# Failures $\left< F \hspace\right>$",
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

        color_index = findfirst(x -> x == freq_bound, filtered_freq_bounds)
        bound_label = ""
        if freq_bound == 0.2
            bound_label = "Narrow bounds"
        elseif freq_bound == 0.3
            bound_label = "Intermediate bounds"
        elseif freq_bound == 1.0
            bound_label = "Wide bounds"
        end

        if (trip_lines == :dynamic &&  trip_nodes == :dynamic)
            if sum_lines_nodes == true
                # scatterlines!(ax_lines_and_nodes, filtered_inertia_values, y_lines + y_nodes, marker = marker, markersize = markersize, label = "f_b=$freq_bound", color = line_colors[color_index])
                # scatterlines!(ax_lines_and_nodes, filtered_inertia_values, y_lines + y_nodes, color = line_colors[color_index], linewidth=4, label = "f_b=$freq_bound")
                scatterlines!(ax_lines_and_nodes, filtered_inertia_values, y_lines + y_nodes, color = line_colors[color_index], linewidth=linewidth, label = "$bound_label")
            else
                scatterlines!(ax_lines_and_nodes, filtered_inertia_values, y_lines, linestyle=:dash, marker = marker, markersize = markersize, label = "f_b=$freq_bound", color = line_colors[color_index])
                band!(ax_lines_and_nodes, filtered_inertia_values, y_lines + err_lines, y_lines - err_lines, transparency=true, color = (line_colors[color_index], opacity))
                scatterlines!(ax_lines_and_nodes, filtered_inertia_values, y_nodes, marker = marker, markersize = markersize, color = line_colors[color_index])
            end
        end
    end
end
M,γ,τ,freq_bound,trip_lines,trip_nodes,init_pert,ensemble_element = RTS_get_network_args_stripped(df_config, 1)

Label(fig_lines_and_nodes[1, 1, TopLeft()], "b", fontsize = 44, font = :bold, padding = (-25, 15, 5, 5))

axislegend(ax_lines_and_nodes, position = :rt, labelsize=labelsize)

# xlims!(ax_lines_and_nodes, 0, 20.5)
# ylims!(ax_lines_and_nodes, 0, 5)

# Save plots
filtered_freq_bounds_str = string(filtered_freq_bounds)
# filtered_freq_bounds_str = "all_frequencies"
# CairoMakie.save(joinpath(MA_DIR, "braessness", "RTS_uebergang_lines+nodes_sumlinesnodes=$sum_lines_nodes,f_b=$filtered_freq_bounds_str,M_left_out=$left_out_inertia_values.pdf"),fig_lines_and_nodes)
# CairoMakie.save(joinpath(MA_DIR, "braessness", "RTS_uebergang_lines+nodes_sumlinesnodes=$sum_lines_nodes,f_b=$filtered_freq_bounds_str,M_left_out=$left_out_inertia_values.png"),fig_lines_and_nodes)
fig_lines_and_nodes