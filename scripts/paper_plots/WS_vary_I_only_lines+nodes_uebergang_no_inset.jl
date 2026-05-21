"""
Watts-Strogatz-Network-Ensemble: Using job array framework. Transition that appears
when varying the frequency bounds. Line and node failures summed.
"""
#  NOTE Check NORMALIZED sum of lines and nodes again (not relevant for figures for paper submission).

include(abspath(@__DIR__, "..", "helpers_jarray.jl"))

using GraphMakie
using Colors, ColorSchemes
using CairoMakie
CairoMakie.activate!()

include(abspath(@__DIR__, "paper_plots_helpers_and_parameters.jl"))

# plotting parameters
create_posprocessing_data = false # set to `false` for fast plotting
sum_lines_nodes = true
normalize = false
# line_colors = [Makie.wong_colors()[3], Makie.wong_colors()[5], Makie.wong_colors()[6]]
line_colors = ["#5C3D99FF", "#56B4E9FF", "#D55E00FF"] # [deep indigo, Hellblau, Orange]
colormap_frequencies = true
opacity = 0.15

# markers
markersize = 15
linewidth = 4
markers_labels = [(:circle, ":circle")]

# exp_name_date = "WS_k=4_exp04_vary_I_only_lines_and_nodes_PIK_HPC_K_=3,N_G=32_20250126_012357.344"
exp_name_date = "WS_k=4_exp04_vary_I_only_lines_and_nodes_PIK_HPC_K_=3,N_G=32_20250321_171511.976"
exp_data_dir = joinpath(RESULTS_DIR, exp_name_date)
left_out_frequencies = [0.005, 0.015, 0.02, 0.025, 0.035, 0.04, 0.045,
    0.05, 0.055, 0.060, 0.065, 0.07, 0.075, 0.08, 0.085, 0.09, 0.095, 0.1,
    0.11, 0.12, 0.13, 0.14, 0.16, 0.17, 0.18, 0.19, 0.2,
    0.21, 0.22, 0.23, 0.24, 0.25, 0.26, 0.27, 0.28, 0.29, 0.3, 0.8]
left_out_inertia_values = []
left_out_β_values = []

################################################################################
###################### Calculate mean and standard error #######################
################################################################################
if create_posprocessing_data == true
    #= # NOTE `postprocess_jarray_data` leads to error for "WS_k=4_exp04_vary_I_only_lines_and_nodes_PIK_HPC_K_=3,N_G=32_20250321_171511.976"
    `preprocess_WS` was changed in order to save power injections. However the data in
    "WS_k=4_exp04_vary_I_only_lines_and_nodes_PIK_HPC_K_=3,N_G=32_20250321_171511.976" were generated before that change.
    =#
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
    title = "Node-line failure model (WS)",
    titlefont = :regular,
    xlabel = L"Inertia $I$ [$s^2$]",
    ylabel = normalize ? "normalized average of failures" : L"# Failures $\left< F \right>$",
)
# hidedecorations!(ax_lines_and_nodes, grid=false)

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
        bound_label = ""
        if freq_bound == 0.01
            bound_label = "Narrow bounds"
        elseif freq_bound == 0.03
            bound_label = "Intermediate bounds"
        elseif freq_bound == 0.15
            bound_label = "Wide bounds"
        end

        if (trip_lines == :dynamic &&  trip_nodes == :dynamic)

            if sum_lines_nodes == true # NOTE For the normalized version this might be wrong s. Schmierzettel S. 23.
                # scatterlines!(ax_lines_and_nodes, filtered_inertia_values, y_lines + y_nodes, marker = marker, markersize = markersize, label = "f_b=$freq_bound,k=$k,β=$β", color = line_colors[color_index])
                scatterlines!(ax_lines_and_nodes, filtered_inertia_values, y_lines + y_nodes, marker = marker, markersize = markersize, linewidth=linewidth, label = "$bound_label", color = line_colors[color_index])
                band!(ax_lines_and_nodes, filtered_inertia_values, y_lines + y_nodes + err_nodes_plus_lines, y_lines + y_nodes - err_nodes_plus_lines, transparency=true, color = (line_colors[color_index], opacity))
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
# lines!(ax_lines_and_nodes, [NaN], [NaN]; label="Damping D=1 [s]", color=:white)

axislegend(ax_lines_and_nodes, position = :rt, labelsize=labelsize)
Label(fig_lines_and_nodes[1, 1, TopLeft()], "a", fontsize = fontsize+labellettersize, font = :bold, padding = (-25, 0, 5, 0))



xlims!(ax_lines_and_nodes, 0, 30.5)
ylims!(ax_lines_and_nodes, 0, 5)


k_str = string(exp_params_dict[:k])
filtered_freq_bounds_str = string(filtered_freq_bounds)
K_str = string(exp_params_dict[:K])

CairoMakie.save(joinpath(MA_DIR, "WS_vary_I_only_uebergang_lines+nodes_sumlinesnodes=$sum_lines_nodes,K=$K_str,k=$k_str,β=$filtered_β_values,f_b=$filtered_freq_bounds_str,M_left_out=$left_out_inertia_values,no_inset.pdf"),fig_lines_and_nodes)
# CairoMakie.save(joinpath(MA_DIR, "WS_vary_I_only_uebergang_lines+nodes_sumlinesnodes=$sum_lines_nodes,K=$K_str,k=$k_str,β=$filtered_β_values,f_b=$filtered_freq_bounds_str,M_left_out=$left_out_inertia_values,no_inset.png"),fig_lines_and_nodes)
fig_lines_and_nodes