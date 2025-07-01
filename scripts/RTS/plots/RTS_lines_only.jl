"""
RTS-GMLC-Testcase: Using job array framework. Plotting lines only.
"""

include(abspath(@__DIR__, "..", "helpers_jarray.jl"))

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
    postprocess_jarray_data(exp_name_date)
end

################################################################################
################################ Plotting  #####################################
################################################################################
inertia_values = exp_params_dict[:inertia_values]
filtered_inertia_values = filter!(x->x ∉ left_out_inertia_values, deepcopy(inertia_values))

fig_lines_only = Figure(size=(800,600),fontsize = fontsize)
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
            scatterlines!(ax_lines_only, filtered_inertia_values, y_lines, label = "", color = Makie.wong_colors()[1], linewidth = 3.5)
        end
    end
end
M,γ,τ,freq_bound,trip_lines,trip_nodes,init_pert,ensemble_element = RTS_get_network_args_stripped(df_config, 1)

# CairoMakie.save(joinpath(MA_DIR, "RTS_lines_only_M_left_out=$left_out_inertia_values.png"),fig_lines_only)
# CairoMakie.save(joinpath(MA_DIR, "RTS_lines_only_M_left_out=$left_out_inertia_values.pdf"),fig_lines_only)
fig_lines_only
