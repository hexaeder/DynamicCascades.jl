include("helpers_jarray.jl")

if ON_YOGA
    using Revise
else # if on PIK-HPC or Pool
    Pkg.instantiate()
    # Pkg.precompile()
end

using GraphMakie
using Colors
using CairoMakie

# exp_name_date = "WS_testrun_params_k=10PIK_HPC_K_=1,N_G=4_20240107_003159.367"
# exp_name_date = "WS_testrun_params_K=6_pool_N_G=2_20240106_021205.818"
# exp_name_date = "WS_testrun_paramsK=9_pool_N_G=2_20240106_020923.414"
# exp_name_date = "WS_testrun_params_K=3_pool_N_G=2_20240106_021759.114"
exp_name_date = "WS_k=4_exp01_PIK_HPC_K_=3,N_G=32_20240128_215811.815"
# exp_name_date = "WS_testrun_params_k=4_narrowPIK_HPC_K_=3,N_G=10_20240119_192800.546"
# exp_name_date = "WS_testrun_params_k=4_widePIK_HPC_K_=3,N_G=10_20240119_192910.398"

# left out frequency values
left_out_frequencies = [0.03]
# left_out_frequencies = [0.001, 0.01, 0.03, 0.05, 0.1, 0.3]
# left_out_frequencies = [0.02, 0.025]

###########
# 0.00100
# 0.0100
# 0.0300
# 0.0500
# 0.100
# 0.300
# 0.500
###########

# left out inertia values
# left_out_inertia_values = [0.01]
left_out_inertia_values = [20.0, 30.0]
# left_out_inertia_values = [ ]

left_out_β_values = []

exp_data_dir = joinpath(RESULTS_DIR, exp_name_date)

################################################################################
###################### Calculate mean and standard error #######################
################################################################################

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
df_avg_error[!, :ensemble_avg_line_failures] .= NaN; df_avg_error[!, :ensemble_avg_node_failures] .= NaN;
df_avg_error[!, :ensemble_SE_line_failures] .= NaN; df_avg_error[!, :ensemble_SE_node_failures] .= NaN;

# DataFrame with failures for ALL ArrayTaskIDs
df_all_failures = deepcopy(df_config)
df_all_failures[!, :norm_avg_line_failures] .= NaN ; df_all_failures[!, :norm_avg_node_failures] .= NaN;

for task_id in df_avg_error.ArrayTaskID
    # loop over all elements of an ensemble
    norm_avg_line_failures_ensemble = Float64[]
    norm_avg_node_failures_ensemble = Float64[]
    for i in 0:num_parameter_combinations:(length(df_config[!,:ArrayTaskID]) - 1)
        # try...catch is for execution of postprocessing while not all jobs have finished
        try
            N,k,β,graph_seed,μ,σ,distr_seed,K,α,M,γ,τ,freq_bound,trip_lines,trip_nodes,init_pert,ensemble_element = get_network_args_stripped(df_config, (task_id + i))

            exp_data = joinpath(RESULTS_DIR, exp_name_date)
            graph_combinations_path = joinpath(exp_data, "k=$k,β=$β")

            failure_mode_string = joinpath(graph_combinations_path, "trip_lines=$trip_lines,trip_nodes=$trip_nodes")
            failure_mode_frequ_bound = joinpath(failure_mode_string, "trip_lines=$trip_lines,trip_nodes=$trip_nodes,freq_bound=$freq_bound")

            filename = string("/", string_network_args(df_config, task_id +i), ".csv")
            df_result = DataFrame(CSV.File(string(failure_mode_frequ_bound, filename)))

            norm_avg_line_failures = df_result[1, :norm_avg_line_failures]
            norm_avg_node_failures = df_result[1, :norm_avg_node_failures]

            push!(norm_avg_line_failures_ensemble, norm_avg_line_failures)
            push!(norm_avg_node_failures_ensemble, norm_avg_node_failures)

            df_all_failures[(task_id + i), :norm_avg_line_failures] = norm_avg_line_failures
            df_all_failures[(task_id + i), :norm_avg_node_failures] = norm_avg_node_failures
        catch
            continue
        end
    end
    # Calculate ensemble_avg and ensemble_standard_error and write to df
    df_avg_error[task_id,:ensemble_avg_line_failures] = mean(norm_avg_line_failures_ensemble)
    df_avg_error[task_id,:ensemble_avg_node_failures] = mean(norm_avg_node_failures_ensemble)
    df_avg_error[task_id,:ensemble_SE_line_failures] = 1 / sqrt(N_ensemble_size) * std(norm_avg_line_failures_ensemble; corrected=true)
    df_avg_error[task_id,:ensemble_SE_node_failures] = 1 / sqrt(N_ensemble_size) * std(norm_avg_node_failures_ensemble; corrected=true)
end


CSV.write(joinpath(RESULTS_DIR, exp_name_date, "avg_error.csv"), df_avg_error)
CSV.write(joinpath(RESULTS_DIR, exp_name_date, "all_failures.csv"), df_all_failures)


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

# Colormaps
# color_map = ColorSchemes.plasma
color_map = ColorSchemes.cividis
# color_map = :blues

# Generate distinct colors based on the filtered_freq_bounds
if length(filtered_freq_bounds) == 1
    line_colors = [RGBA{Float64}(1.0,0.9169,0.2731,1.0)]
else
    line_colors = distinct_colors(color_map, filtered_freq_bounds)
end
# markers
markers_labels = [
    (:circle, ":circle"),
    (:rect, ":rect"),
    (:utriangle, ":utriangle"),
]

function create_figs(failure_modes)
    fig_lines_only = Figure(); fig_nodes_only = Figure(); fig_lines_and_nodes= Figure();
    ax_lines_only = Axis(fig_lines_only[1, 1]); ax_nodes_only = Axis(fig_nodes_only[1, 1]); ax_lines_and_nodes = Axis(fig_lines_and_nodes[1, 1])
    for i in exp_params_dict[:failure_modes]
        if i == [:dynamic, :none]
            fig_lines_only = Figure(fontsize = 28)
            ax_lines_only = Axis(fig_lines_only[1, 1],
                title = "Line failures",
                # titlesize = 30,
                xlabel = "scaling factor of inertia",
                # xlabelsize = 30,
                ylabel = "normalized average of line failures",
                # ylabelsize = 30
            )
        elseif i == [:none, :dynamic]
            fig_nodes_only = Figure(fontsize = 28)
            ax_nodes_only = Axis(fig_nodes_only[1, 1],
                title = "Node failures",
                xlabel = "scaling factor of inertia",
                ylabel = "normalized average of node failures",
            )
        elseif i == [:dynamic, :dynamic]
            fig_lines_and_nodes = Figure(fontsize = 28)
            ax_lines_and_nodes = Axis(fig_lines_and_nodes[1, 1],
                title = "Line and node failures",
                xlabel = "scaling factor of inertia",
                ylabel = "normalized average of failures",
            )
        end
    end
    return fig_lines_only, ax_lines_only, fig_nodes_only, ax_nodes_only, fig_lines_and_nodes, ax_lines_and_nodes
end

# Create figures depending on the modes (loop).
failure_modes = exp_params_dict[:failure_modes]
β_vals = exp_params_dict[:β]
fig_lines_only, ax_lines_only, fig_nodes_only, ax_nodes_only, fig_lines_and_nodes, ax_lines_and_nodes = create_figs(failure_modes)

df_avg_error = DataFrame(CSV.File(joinpath(RESULTS_DIR, exp_name_date, "avg_error.csv")))
inertia_values = exp_params_dict[:inertia_values]


#= Different inertia values for: Ensemble average over normalized average of
failures (the latter for a single network) =#
y_lines = Float64[]; y_nodes = Float64[]
#= Different inertia values for: Ensemble standard error over normalized average
of failures (the latter for a single network) =#
err_lines = Float64[]; err_nodes = Float64[]
for task_id in df_avg_error.ArrayTaskID # TODO renane variables: this is not an ArrayTaskID in the strict sense but an average over task IDs
    #= Empty arrays (after all inertia values of one configuration pushed to array)
    The entries in df_avg_error are ordered accordningly.=#
    if ((task_id-1) % length(inertia_values)) == 0
        y_lines = Float64[]; y_nodes = Float64[]
        err_lines = Float64[]; err_nodes = Float64[]
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
    push!(y_lines, df_avg_error[task_id, :ensemble_avg_line_failures])
    push!(y_nodes, df_avg_error[task_id, :ensemble_avg_node_failures])
    push!(err_lines, df_avg_error[task_id, :ensemble_SE_line_failures])
    push!(err_nodes, df_avg_error[task_id, :ensemble_SE_node_failures])

    # Only plot if all inertia values are pushed to `y_lines`, `y_nodes`, `err_lines`, `err_nodes`
    if M == maximum(filtered_inertia_values)
        # frequency argument first for a nice order in the legend
        N,k,β,graph_seed,μ,σ,distr_seed,K,α,M,γ,τ,freq_bound,trip_lines,trip_nodes,init_pert,ensemble_element = get_network_args_stripped(df_config, task_id)

        marker_index = findfirst(x -> x == β, β_vals)
        marker = markers_labels[marker_index][1]
        marker_label = markers_labels[marker_index][2]
        color_index = findfirst(x -> x == freq_bound, filtered_freq_bounds)
        if (trip_lines == :dynamic &&  trip_nodes == :none)
            if freq_bound == filtered_freq_bounds[1]
                scatterlines!(ax_lines_only, filtered_inertia_values, y_lines, marker = marker, label = "k=$k,β=$β")
                errorbars!(ax_lines_only, filtered_inertia_values, y_lines, err_lines, color = :black,  whiskerwidth = 10)
            end
        elseif (trip_lines == :none &&  trip_nodes == :dynamic)
            scatterlines!(ax_nodes_only, filtered_inertia_values, y_nodes, marker = marker, label = "f_b=$freq_bound,k=$k,β=$β", color = line_colors[color_index])
            errorbars!(ax_nodes_only, filtered_inertia_values, y_nodes, err_nodes, color = :black,  whiskerwidth = 10)
        elseif (trip_lines == :dynamic &&  trip_nodes == :dynamic)
            scatterlines!(ax_lines_and_nodes, filtered_inertia_values, y_lines, linestyle=:dash, marker = marker, label = "f_b=$freq_bound,k=$k,β=$β", color = line_colors[color_index])
            errorbars!(ax_lines_and_nodes, filtered_inertia_values, y_lines, err_lines, color = :black,  whiskerwidth = 10)
            scatterlines!(ax_lines_and_nodes, filtered_inertia_values, y_nodes, marker = marker, color = line_colors[color_index])
            errorbars!(ax_lines_and_nodes, filtered_inertia_values, y_nodes, err_nodes, color = :black,  whiskerwidth = 10)
        end
    end
end
N,k,β,graph_seed,μ,σ,distr_seed,K,α,M,γ,τ,freq_bound,trip_lines,trip_nodes,init_pert,ensemble_element = get_network_args_stripped(df_config, 1)
#  See https://docs.makie.org/stable/reference/blocks/legend/
lines!(ax_lines_and_nodes, [NaN], [NaN]; label="node failures: ___ (solid) \nline failures:   ----- (dashed)", color=:white, linewidth=3)
axislegend(ax_lines_only, position = :rb, labelsize=10)
axislegend(ax_nodes_only, position = :rt, labelsize=10)
axislegend(ax_lines_and_nodes, position = :rt, labelsize=10)

# text!(ax_lines_and_nodes, 5.0, 0.02, text = "node failures: solid lines ___ \n line failures: dashed lines -----", align = (:center, :center), textsize=25)

# fig_lines_only
# fig_nodes_only
# fig_lines_and_nodes

# Save plots
k_str = string(exp_params_dict[:k])
filtered_freq_bounds_str = string(filtered_freq_bounds)

K_str = string(exp_params_dict[:K])

CairoMakie.save(joinpath(exp_data_dir, "lines_only_K=$K_str,k=$k_str,β=$filtered_β_values,f_b=$filtered_freq_bounds_str,M_left_out=$left_out_inertia_values.pdf"),fig_lines_only)
CairoMakie.save(joinpath(exp_data_dir, "nodes_only_K=$K_str,k=$k_str,β=$filtered_β_values,f_b=$filtered_freq_bounds_str,M_left_out=$left_out_inertia_values.pdf"),fig_nodes_only)
CairoMakie.save(joinpath(exp_data_dir, "lines+nodes_K=$K_str,k=$k_str,β=$filtered_β_values,f_b=$filtered_freq_bounds_str,M_left_out=$left_out_inertia_values.pdf"),fig_lines_and_nodes)

# # NOTE Further plotting options:
# # Placing the legend besides the coordinate system.
# leg_lines = Legend(fig_lines_only[1, 2], ax_lines_only, "Headline Legend", framevisible = false)
# leg_nodes = Legend(fig_nodes_only[1, 2], ax_nodes_only)
# leg_lines_and_nodes = Legend(fig_lines_and_nodes[1, 2], ax_lines_and_nodes)
# axislegend("Legend headline", position = :rt)
# # Template on how to put all plots in one single figure:
# # https://docs.makie.org/stable/tutorials/layout-tutorial/#panel_a
