"""
Watts-Strogatz-Network-Ensemble: Using job array framework. Transition that appears
when varying the frequency bounds. Line and node failures summed.
"""
#  NOTE NOTE NOTE NOTE NOTE NOTE NOTE NOTE NOTE NOTE NOTE NOTE Check NORMALIZED sum of lines and nodes again.


include(abspath(@__DIR__, "..", "helpers_jarray.jl"))

using GraphMakie
using Colors, ColorSchemes
using CairoMakie


# plotting parameters
create_posprocessing_data = false # set to `false` for fast plotting
sum_lines_nodes = true
normalize = false
line_colors = [Makie.wong_colors()[3], Makie.wong_colors()[5], Makie.wong_colors()[6]]
colormap_frequencies = true
opacity = 0.15
fontsize = labelsize = 26
# markers
markersize = 15
markers_labels = [(:circle, ":circle")]

exp_name_date = "WS_k=4_exp06_3_I_over_D_lines_and_nodes_PIK_HPC_K_=3,N_G=32_20250125_134157.625"
exp_data_dir = joinpath(RESULTS_DIR, exp_name_date)
left_out_frequencies = [0.005, 0.01, 0.015, 0.025, 0.03, 0.035, 0.14, 0.15, 0.16, 0.3, 0.5, 0.8]
left_out_frequencies = [0.025, 0.03, 0.035, 0.14, 0.15, 0.16, 0.3, 0.5, 0.8]
left_out_inertia_values = []
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
################################ Plotting  #####################################
################################################################################
inertia_values = exp_params_dict[:inertia_values]
freq_bounds = exp_params_dict[:freq_bounds]

filtered_freq_bounds = filter!(x->x ∉ left_out_frequencies, deepcopy(freq_bounds))
filtered_inertia_values = filter!(x->x ∉ left_out_inertia_values, deepcopy(inertia_values))
filtered_β_values = filter!(x->x ∉ left_out_β_values, deepcopy(exp_params_dict[:β]))

fig_lines_and_nodes = Figure(size=(800,600),fontsize = fontsize)
ax_lines_and_nodes = Axis(fig_lines_and_nodes[1, 1],
    # title = sum_lines_nodes ? "Summed line and node failures" : "Line and node failures",
    xlabel = L"Inertia I [$s^2$]",
    ylabel = normalize ? "normalized average of failures" : L"Averaged failures $N_{fail}$",
)
# hidedecorations!(ax_lines_and_nodes, grid=false)

# Create inset
# The inset axis
inset_ax = Axis(fig_lines_and_nodes[1, 1],
    # xlabel = L"Inertia I [$s^2$]",
    # ylabel = normalize ? "normalized average of failures" : L"Averaged failures $N_{fail}$",
    width=Relative(0.5),
    height=Relative(0.5),
    halign=0.95,
    valign=0.95,
    backgroundcolor=(:grey, 0.05),
    xgridvisible = false,
    ygridvisible = false,
    )


# hidedecorations!(inset_ax)
xlims!(inset_ax, 0, 30)
ylims!(inset_ax, 0, 1)


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

        marker_index = findfirst(x -> x == β, β_vals)
        marker = markers_labels[marker_index][1]
        marker_label = markers_labels[marker_index][2]
        color_index = colormap_frequencies ? findfirst(x -> x == freq_bound, filtered_freq_bounds) : marker_index

        if (trip_lines == :dynamic &&  trip_nodes == :dynamic)
            if sum_lines_nodes == true # NOTE For the normalized version this might be wrong s. Schmierzettel S. 23.
                scatterlines!(ax_lines_and_nodes, filtered_inertia_values, y_lines + y_nodes, marker = marker, markersize = markersize, label = "f_b=$freq_bound,k=$k,β=$β", color = line_colors[color_index])
                band!(ax_lines_and_nodes, filtered_inertia_values, y_lines + y_nodes + err_nodes_plus_lines, y_lines + y_nodes - err_nodes_plus_lines, transparency=true, color = (line_colors[color_index], opacity))
                # inset plot
                scatterlines!(inset_ax, filtered_inertia_values, y_lines + y_nodes, marker = marker, markersize = markersize, color = line_colors[color_index])
                band!(inset_ax, filtered_inertia_values, y_lines + y_nodes + err_nodes_plus_lines, y_lines + y_nodes - err_nodes_plus_lines, transparency=true, color = (line_colors[color_index], opacity))
            else
                scatterlines!(ax_lines_and_nodes, filtered_inertia_values, y_lines, linestyle=:dash, marker = marker, markersize = markersize, label = "f_b=$freq_bound,k=$k,β=$β", color = line_colors[color_index])
                band!(ax_lines_and_nodes, filtered_inertia_values, y_lines + err_lines, y_lines - err_lines, transparency=true, color = (line_colors[color_index], opacity))
                scatterlines!(ax_lines_and_nodes, filtered_inertia_values, y_nodes, marker = marker, markersize = markersize, color = line_colors[color_index])
                band!(ax_lines_and_nodes, filtered_inertia_values, y_nodes + err_nodes, y_nodes - err_nodes, transparency=true, color = (line_colors[color_index], opacity))
            end
        end
    end
end
N,k,β,graph_seed,μ,σ,distr_seed,K,α,M,γ,τ,freq_bound,trip_lines,trip_nodes,init_pert,ensemble_element = get_network_args_stripped(df_config, 1)
lines!(ax_lines_and_nodes, [NaN], [NaN]; label=L"Damping $D=I$ [$s$]", color=:white)
axislegend(ax_lines_and_nodes, position = (0.07,0.99), labelsize=labelsize-5)

k_str = string(exp_params_dict[:k])
filtered_freq_bounds_str = string(filtered_freq_bounds)
K_str = string(exp_params_dict[:K])

CairoMakie.save(joinpath(MA_DIR, "WS", "WS_I_over_D_uebergang_lines+nodes_sumlinesnodes=$sum_lines_nodes,K=$K_str,k=$k_str,β=$filtered_β_values,f_b=$filtered_freq_bounds_str,M_left_out=$left_out_inertia_values.png"),fig_lines_and_nodes)
CairoMakie.save(joinpath(MA_DIR, "WS", "WS_I_over_D_uebergang_lines+nodes_sumlinesnodes=$sum_lines_nodes,K=$K_str,k=$k_str,β=$filtered_β_values,f_b=$filtered_freq_bounds_str,M_left_out=$left_out_inertia_values.pdf"),fig_lines_and_nodes)
fig_lines_and_nodes
