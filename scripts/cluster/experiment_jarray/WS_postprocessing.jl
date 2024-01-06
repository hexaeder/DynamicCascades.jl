include("helpers_jarray.jl")

if ON_YOGA
    using Revise
else # if on PIK-HPC or Pool
    Pkg.instantiate()
    # Pkg.precompile()
end

using LinearAlgebra
print("Number of threads before setting"); print(LinearAlgebra.BLAS.get_num_threads()); print("\n")
BLAS.set_num_threads(1)
print("Number of threads after setting"); print(LinearAlgebra.BLAS.get_num_threads()); print("\n")

using DynamicCascades
using NetworkDynamics
using Graphs
using MetaGraphs
using Unitful
using Statistics
using Dates
using DataFrames
using CSV
using Serialization
using GraphMakie
using Colors
using CairoMakie

exp_name_date = "WS_testrun_plots_N_G=3_20240103_202725.758"
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

# DataFrame with failures for all ArrayTaskIDs
df_all_failures = deepcopy(df_config)
df_all_failures[!, :norm_avg_line_failures] .= NaN ; df_all_failures[!, :norm_avg_node_failures] .= NaN;

for task_id in df_avg_error.ArrayTaskID
    # loop over all elements of an ensemble
    norm_avg_line_failures_ensemble = Float64[]
    norm_avg_node_failures_ensemble = Float64[]
    for i in 0:num_parameter_combinations:(length(df_config[!,:ArrayTaskID]) - 1)
        try
            N,k,β,graph_seed,μ,σ,distr_seed,K,α,M,γ,τ,freq_bound,trip_lines,trip_nodes,init_pert,ensemble_element = get_network_args_stripped(df_config, (task_id + i))

            exp_data = joinpath(RESULTS_DIR, exp_name_date)
            graph_combinations_path = joinpath(exp_data, "k=$k,β=$β")

            failure_mode_string = joinpath(graph_combinations_path, "trip_lines=$trip_lines,trip_nodes=$trip_nodes")
            failure_mode_frequ_bound = joinpath(failure_mode_string, "trip_lines=$trip_lines,trip_nodes=$trip_nodes,freq_bound=$freq_bound")

            filename = string("/", string_network_args(df_config, task_id), ".csv")
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
freq_bounds = exp_params_dict[:freq_bounds]
fig_lines_only, ax_lines_only, fig_nodes_only, ax_nodes_only, fig_lines_and_nodes, ax_lines_and_nodes = create_figs(failure_modes)

df_avg_error = DataFrame(CSV.File(joinpath(RESULTS_DIR, exp_name_date, "avg_error.csv")))
inertia_values = exp_params_dict[:inertia_values]
# lines!([NaN], [NaN]; label="Bounds", color=:white, linewidth=3) # headline legend
for task_id in df_avg_error.ArrayTaskID
    #= Empty arrays (after all inertia values of one configuration pushed to array)
    The entries in df_avg_error are ordered accordningly.=#
    if ((task_id-1) % length(inertia_values)) == 0
        #= Different inertia values for: Ensemble average over normalized average of
        failures (the latter for a single network) =#
        y_lines = Float64[]; y_nodes = Float64[]
        #= Different inertia values for: Ensemble standard error over normalized average
        of failures (the latter for a single network) =#
        err_lines = Float64[]; err_nodes = Float64[]
    end

    # Read out ensemble_avg and ensemble_standard_error
    push!(y_lines, df_avg_error[task_id, :ensemble_avg_line_failures])
    push!(y_nodes, df_avg_error[task_id, :ensemble_avg_node_failures])
    push!(err_lines, df_avg_error[task_id, :ensemble_SE_line_failures])
    push!(err_nodes, df_avg_error[task_id, :ensemble_SE_node_failures])

    if (task_id % length(inertia_values)) == 0
        # frequency argument first for a nice order in the legend
        N,k,β,graph_seed,μ,σ,distr_seed,K,α,M,γ,τ,freq_bound,trip_lines,trip_nodes,init_pert,ensemble_element = get_network_args_stripped(df_config, task_id)

        if (trip_lines == :dynamic &&  trip_nodes == :none)
            if (task_id % (length(inertia_values) * length(freq_bounds))) == 0
                scatterlines!(ax_lines_only, inertia_values, y_lines, label = "k=$k,β=$β")
                errorbars!(ax_lines_only, inertia_values, y_lines, err_lines, color = :black,  whiskerwidth = 10)
            end
        elseif (trip_lines == :none &&  trip_nodes == :dynamic)
            scatterlines!(ax_nodes_only, inertia_values, y_nodes, label = "f_b=$freq_bound,k=$k,β=$β")
            errorbars!(ax_nodes_only, inertia_values, y_nodes, err_nodes, color = :black,  whiskerwidth = 10)
        elseif (trip_lines == :dynamic &&  trip_nodes == :dynamic)
            scatterlines!(ax_lines_and_nodes, inertia_values, y_lines, label = "f_b=$freq_bound,k=$k,β=$β")
            errorbars!(ax_lines_and_nodes, inertia_values, y_lines, err_lines, color = :black,  whiskerwidth = 10)
            scatterlines!(ax_lines_and_nodes, inertia_values, y_nodes, linestyle=:dash)
            errorbars!(ax_lines_and_nodes, inertia_values, y_nodes, err_nodes, color = :black,  whiskerwidth = 10)
        end
    end
end

#  See https://docs.makie.org/stable/reference/blocks/legend/
axislegend(ax_lines_only, position = :rt)
axislegend(ax_nodes_only, position = :rt)
axislegend(ax_lines_and_nodes, position = :rt)
text!(ax_lines_and_nodes, 0.5, 0.02, text = "node failures: solid lines ___ \n line failures: dashed lines -----", align = (:center, :center), textsize=25)

# fig_lines_only
# fig_nodes_only
# fig_lines_and_nodes

# Save plots
k_str = string(exp_params_dict[:k])
β_str = string(exp_params_dict[:β])
freq_bounds_str = string(exp_params_dict[:freq_bounds])

CairoMakie.save(joinpath(exp_data_dir, "lines_only_k=$k_str,β=$β_str,f_b=$freq_bounds_str.pdf"),fig_lines_only)
CairoMakie.save(joinpath(exp_data_dir, "nodes_only_k=$k_str,β=$β_str,f_b=$freq_bounds_str.pdf"),fig_nodes_only)
CairoMakie.save(joinpath(exp_data_dir, "lines+nodes_k=$k_str,β=$β_str,f_b=$freq_bounds_str.pdf"),fig_lines_and_nodes)


# # NOTE Further plotting options:
# # Placing the legend besides the coordinate system.
# leg_lines = Legend(fig_lines_only[1, 2], ax_lines_only, "Headline Legend", framevisible = false)
# leg_nodes = Legend(fig_nodes_only[1, 2], ax_nodes_only)
# leg_lines_and_nodes = Legend(fig_lines_and_nodes[1, 2], ax_lines_and_nodes)
# axislegend("Legend headline", position = :rt)
# # Template on how to put all plots in one single figure:
# # https://docs.makie.org/stable/tutorials/layout-tutorial/#panel_a
#
# # Color coding of lines
# using Colors, ColorSchemes
# using Printf
#
# ang_freq_bounds = [0.24, 0.32, 0.40, 0.48, 0.64, 0.80, 0.95, 1.27, 1.59, 1.91]
#
# norm_values = (ang_freq_bounds .- minimum(ang_freq_bounds)) ./ (maximum(ang_freq_bounds) - minimum(ang_freq_bounds))
#
# # Function to generate distinct colors based on a color map and normalization
# function distinct_colors(color_map, values)
#     norm_values = (values .- minimum(values)) ./ (maximum(values) - minimum(values))
#
#     colors = [cgrad(color_map, 101; categorical = true, rev=true)[Int(ceil(i*100)+1)] for i in norm_values]
#     return colors
# end
# # color_map = ColorSchemes.plasma
# color_map = ColorSchemes.cividis
# # color_map = :blues
#
# # Generate distinct colors based on the ang_freq_bounds
# line_colors = distinct_colors(color_map, ang_freq_bounds)
#
# for (i, dir, ang_freq) in zip(1:length(line_colors),directories, ang_freq_bounds)
#     directory = string(MA_DATA, dir)
#
#     x = df_inertia_vs_failures_nodes.scale_inertia_values
#     y_nodes = df_inertia_vs_failures_nodes.rel_failures
#     y_lines = df_inertia_vs_failures_lines.rel_failures
#
#     scatterlines!(x, y_nodes, label = "$ang_freq", color = line_colors[i])
#     scatterlines!(x, y_lines, color = line_colors[i], linestyle=:dash )
# end
