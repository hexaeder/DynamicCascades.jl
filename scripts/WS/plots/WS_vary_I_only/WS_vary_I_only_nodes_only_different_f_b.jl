"""
Watts-Strogatz-Network-Ensemble: Using job array framework. Only varying frequency
bound f_b and only considering node failures.
"""

include(abspath(@__DIR__, "..", "helpers_jarray.jl"))

using GraphMakie
using Colors
using CairoMakie


# plotting parameters
create_posprocessing_data = true # set to `false` for fast plotting
sum_lines_nodes = true
normalize = false
custom_colors = true
predefined_colors = [Makie.wong_colors()[1], Makie.wong_colors()[2], Makie.wong_colors()[3], Makie.wong_colors()[4]]  # https://docs.makie.org/stable/explanations/colors/
colormap_frequencies = true
opacity = 0.3
fontsize = labelsize = 26
# markers
markersize = 15
markers_labels = [
    (:utriangle, ":utriangle"),
    (:rect, ":rect"),
    (:star5, "star5"),
    (:circle, ":circle"),
]

exp_name_date = "WS_k=4_exp03_2_vary_I_only_nodes_PIK_HPC_K_=3,N_G=32_20250124_123411.73"
exp_data_dir = joinpath(RESULTS_DIR, exp_name_date)
left_out_frequencies = [0.010, 0.015, 0.020, 0.025, 0.035, 0.040,
    0.045, 0.050, 0.055, 0.060, 0.065, 0.070, 0.075, 0.080, 0.085, 0.090, 0.095, 0.100,
    0.110, 0.120, 0.130, 0.140, 0.150, 0.160, 0.170, 0.180, 0.190, 0.200,
    0.210, 0.220, 0.230, 0.240, 0.250, 0.260, 0.270, 0.280, 0.290, 0.300, 0.800]
left_out_inertia_values = []
left_out_β_values = [0.25, 0.5]

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

# Color coding of lines
using Colors, ColorSchemes

norm_values = (filtered_freq_bounds .- minimum(filtered_freq_bounds)) ./ (maximum(filtered_freq_bounds) - minimum(filtered_freq_bounds))

# Function to generate distinct colors based on a color map and normalization
function distinct_colors(color_map, values)
    norm_values = (values .- minimum(values)) ./ (maximum(values) - minimum(values))

    colors = [cgrad(color_map, 101; categorical = true, rev=true)[Int(ceil(i*100)+1)] for i in norm_values]
    return colors
end

color_map =  [Makie.wong_colors()[1], Makie.wong_colors()[5]]


# Generate distinct colors based on the filtered_freq_bounds
if colormap_frequencies
    line_colors = distinct_colors(color_map, filtered_freq_bounds)
    if length(filtered_freq_bounds) == 1
        # line_colors = [RGBA{Float64}(0.0,0.1262,0.3015,1.0)]
        line_colors = [Makie.wong_colors()[1]]
    end
elseif custom_colors
    line_colors = predefined_colors
end

fig_nodes_only = Figure(size=(800,600),fontsize = fontsize)
ax_nodes_only = Axis(fig_nodes_only[1, 1],
    title = "",
    xlabel = L"Inertia I [$s^2$]",
    ylabel = normalize ? "normalized average of node failures" : L"Averaged node failures $N_{fail}^N$",
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

        # if colormap_frequencies
        #     color_index = findfirst(x -> x == freq_bound, filtered_freq_bounds)
        # end

        if (trip_lines == :none &&  trip_nodes == :dynamic)
            scatterlines!(ax_nodes_only, filtered_inertia_values, y_nodes, marker = marker,  markersize = markersize, label = "f_b=$freq_bound,k=$k,β=$β", color = line_colors[color_index])
            band!(ax_nodes_only, filtered_inertia_values, y_nodes + err_nodes, y_nodes - err_nodes, transparency=true, color = (line_colors[color_index], opacity))
        end
    end
end
N,k,β,graph_seed,μ,σ,distr_seed,K,α,M,γ,τ,freq_bound,trip_lines,trip_nodes,init_pert,ensemble_element = get_network_args_stripped(df_config, 1)
lines!(ax_nodes_only, [NaN], [NaN]; label="Damping D=1 [s]", color=:white)
axislegend(ax_nodes_only, position = :rt, labelsize=labelsize)

k_str = string(exp_params_dict[:k])
filtered_freq_bounds_str = string(filtered_freq_bounds)
K_str = string(exp_params_dict[:K])

CairoMakie.save(joinpath(MA_DIR, "WS", "WS_vary_I_only_nodes_only_K=$K_str,k=$k_str,β=$filtered_β_values,f_b=$filtered_freq_bounds_str,M_left_out=$left_out_inertia_values.png"),fig_nodes_only)
CairoMakie.save(joinpath(MA_DIR, "WS", "WS_vary_I_only_nodes_only_K=$K_str,k=$k_str,β=$filtered_β_values,f_b=$filtered_freq_bounds_str,M_left_out=$left_out_inertia_values.pdf"),fig_nodes_only)
fig_nodes_only
