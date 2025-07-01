"""
Watts-Strogatz-Network-Ensemble: Using job array framework.
Heatmap: x: frequency bound f_b, y: optimal Inertia I
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
heatmap_logscale = true
opacity = 0.3
fontsize = labelsize = 24
# markers
markersize = 15

exp_name_date = "WS_k=4_exp04_vary_I_only_lines_and_nodes_PIK_HPC_K_=3,N_G=32_20250126_012357.344"
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

# Save plots
k_str = string(exp_params_dict[:k])
filtered_freq_bounds_str = string(filtered_freq_bounds)
K_str = string(exp_params_dict[:K])

if length(filtered_freq_bounds) > 1
    # create minimal failures (optimal inertia) vs. frequency bound f_b ########
    fig_opt_inertia_vs_fb = Figure(size=(800,600),fontsize = (fontsize-3))
    ax_opt_inertia_vs_fb = Axis(fig_opt_inertia_vs_fb[1, 1], xticklabelrotation=π/2,
        title = "",
        xlabel = "Frequency bound f_b [Hz]", # frequency bound f_b
        xlabelsize = (fontsize + 5),
        ylabel = L"Inertia $I_{min}$ [$s^2$]", # inertia value associated with minimum of failures
        ylabelsize = (fontsize + 5),
    )

    ax_opt_inertia_vs_fb.xticks = filtered_freq_bounds
    ax_opt_inertia_vs_fb.yticks = [0.2, 3.0, 5.0, 10.0, 20.0, 30.0]

    for i in 1:length(filtered_freq_bounds)
        scatter!(ax_opt_inertia_vs_fb, filtered_freq_bounds[i], opt_inertia[i], color = Makie.wong_colors()[1], markersize = markersize)
    end

    # vspan!([0.0025, 0.0125, 0.0775], [0.0125, 0.0775, 0.155], color = [(:grey, c) for c in [0.1, 0.3, 0.5]])
    vspan!([0.0025, 0.0125, 0.0975], [0.0125, 0.0975, 0.155], color = [(:grey, c) for c in [0.1, 0.3, 0.5]])
    # text!(ax_opt_inertia_vs_fb, 0.013, 8., text = "I", align = (:center, :center), rotation=π/2, textsize=35, color=:pink)
    text!(ax_opt_inertia_vs_fb, 0.007, 25., text = L"I", align = (:center, :center), textsize=45)
    text!(ax_opt_inertia_vs_fb, 0.035, 25., text = L"II", align = (:center, :center), textsize=45)
    text!(ax_opt_inertia_vs_fb, 0.11, 25., text = L"III", align = (:center, :center), textsize=45)

end
CairoMakie.save(joinpath(MA_DIR, "WS", "WS_vary_I_only_optimal_inertia_sumlinesnodes=$sum_lines_nodes,K=$K_str,k=$k_str,β=$filtered_β_values,f_b_left_out=[_],M_left_out=$left_out_inertia_values.pdf"),fig_opt_inertia_vs_fb)
CairoMakie.save(joinpath(MA_DIR, "WS", "WS_vary_I_only_optimal_inertia_sumlinesnodes=$sum_lines_nodes,K=$K_str,k=$k_str,β=$filtered_β_values,f_b_left_out=[_],M_left_out=$left_out_inertia_values.png"),fig_opt_inertia_vs_fb)
fig_opt_inertia_vs_fb
