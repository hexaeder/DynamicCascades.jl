"""
Watts-Strogatz-Network-Ensemble: Using job array framework.
Heatmap: x: frequency bound f_b, y: Inertia I, z: log(sum of line + node failres + 1)
"""
#  NOTE NOTE NOTE NOTE NOTE NOTE NOTE NOTE NOTE NOTE NOTE NOTE Check normalized sum of lines and nodes again.

include(abspath(@__DIR__, "..", "..", "..", "helpers_jarray.jl"))

using GraphMakie
using Colors
using CairoMakie
CairoMakie.activate!()

# plotting parameters
create_posprocessing_data = false # set to `false` for fast plotting
sum_lines_nodes = true
normalize = false
heatmap_logscale = true
opacity = 0.3
fontsize = labelsize = 24
# markers
markersize = 15

exp_name_date = "WS_k=4_exp04_vary_I_only_lines_and_nodes_PIK_HPC_K_=3,N_G=32_20250321_171511.976"
exp_data_dir = joinpath(RESULTS_DIR, exp_name_date)
left_out_frequencies = [
    0.16, 0.17, 0.18, 0.19, 0.2, 0.21, 0.22, 0.23, 0.24, 0.25, 0.26, 0.27, 0.28, 0.29, 0.3, 0.8]
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
    fig_hm = Figure(size=(800,600),fontsize = (fontsize-3))
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
    lines!(ax_hm, [NaN], [NaN]; label=L"Damping $D=1$ [$s$]", color=:white)
    axislegend(ax_hm, position = :rt, labelsize=(labelsize-8))
    Colorbar(fig_hm[:, end+1], hm, label = normalize ? L"normalized $N_{fail}$" : (heatmap_logscale ? L"$\log(N_{fail}+1)$" : L"$N_{fail}$"))

    ax_hm.xticks = filtered_freq_bounds
    ax_hm.xlabelpadding = 15
    ax_hm.yticks = [1.0, 3.0, 5.0, 7.5, 10.0, 20.0, 30.0]
end
# CairoMakie.save(joinpath(MA_DIR, "WS", "WS_vary_I_only_heatmap_log=$heatmap_logscale,sumlinesnodes=$sum_lines_nodes,K=$K_str,k=$k_str,β=$filtered_β_values.pdf"),fig_hm)
# CairoMakie.save(joinpath(MA_DIR, "WS", "WS_vary_I_only_heatmap_log=$heatmap_logscale,sumlinesnodes=$sum_lines_nodes,K=$K_str,k=$k_str,β=$filtered_β_values.png"),fig_hm)
fig_hm
