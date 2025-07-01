"""
RTS-GMLC-Testcase: Not using job array framework. Plotting nodes only.
"""

# NOTE Many packages below not needed.
using DynamicCascades
using NetworkDynamics
using Graphs
using MetaGraphs
using Unitful
using Statistics
using GraphMakie
using Colors
using DynamicCascades: PLOT_DIR
using CairoMakie
using Dates
using DataFrames
using CSV
using Colors, ColorSchemes
using Printf

show_title = false
fontsize = labelsize = 28
line_colors = [Makie.wong_colors()[1], Makie.wong_colors()[3], Makie.wong_colors()[5], Makie.wong_colors()[6]]

############################### node failures only #############################
directories = ["/results_plotting/node_failures_only/20230729_034424.446inertia_vs_line+node_failures_f_bound=0.24",
            "/results_plotting/node_failures_only/20230729_034425.247inertia_vs_line+node_failures_f_bound=0.32",
            # "/results_plotting/node_failures_only/20230729_034431.966inertia_vs_line+node_failures_f_bound=0.4",
            # "/results_plotting/node_failures_only/20230729_034435.458inertia_vs_line+node_failures_f_bound=0.48",
            "/results_plotting/node_failures_only/20230730_120321.717inertia_vs_line+node_failures_f_bound=0.64",
            # "/results_plotting/node_failures_only/20230730_121128.261inertia_vs_line+node_failures_f_bound=0.8",
            # "/results_plotting/node_failures_only/20230730_121128.342inertia_vs_line+node_failures_f_bound=0.95",
            # "/results_plotting/node_failures_only/20230730_121128.333inertia_vs_line+node_failures_f_bound=1.27",
            "/results_plotting/node_failures_only/20230730_121128.315inertia_vs_line+node_failures_f_bound=1.59"
            # "/results_plotting/node_failures_only/20230730_121159.099inertia_vs_line+node_failures_f_bound=1.91"
            ]

ang_freq_bounds = [0.24, 0.32, 0.64, 1.59]

# plot data
fig = Figure(size=(800,600),fontsize = fontsize)
Axis(fig[1, 1],
    title = show_title ? "Variation of frequency bounds" : "",
    xlabel = "Scaling factor of inertia I",
    ylabel = L"Averaged node failures $N_{fail}^N$",
)

network = import_system(:rtsgmlc; damping = 0.1u"s", tconst = 0.01u"s")
nd, = nd_model(network)
ω_state_idxs = idx_containing(nd, "ω")
gen_node_idxs = map(s -> parse(Int, String(s)[4:end]), nd.syms[ω_state_idxs])
# lines!([NaN], [NaN]; label="Bounds", color=:white, linewidth=3)

for (i, dir, ang_freq) in zip(1:length(line_colors),directories, ang_freq_bounds)
    directory = string(MA_DATA, dir)
    df_all_failures_nodes = DataFrame(CSV.File(string(directory,"/all_failures_nodes.csv")))

    # calculate relative number of failures
    rel_number_failures = Float64[] # NOTE this is actually the mean
    scale_inertia_values = names(df_all_failures_nodes)
    for scale_inertia in scale_inertia_values
        number_failures_nodes = df_all_failures_nodes[!, scale_inertia]
        # NOTE TODO `mean(number_failures_nodes)` is not correct as long as line 11 is not included
        # # normalized
        # push!(rel_number_failures, mean(number_failures_nodes)/length(gen_node_idxs))
        # not normalized
        push!(rel_number_failures, mean(number_failures_nodes))
    end
    df_inertia_vs_failures_nodes = DataFrame("scale_inertia_values" => scale_inertia_values, "rel_failures" => rel_number_failures)
    CSV.write(string(directory,"/inertia_vs_failures_nodes.csv"), df_inertia_vs_failures_nodes)

    # load data
    df_inertia_vs_failures_nodes = DataFrame(CSV.File(string(directory,"/inertia_vs_failures_nodes.csv")))

    x = df_inertia_vs_failures_nodes.scale_inertia_values
    y_nodes = df_inertia_vs_failures_nodes.rel_failures

    scatterlines!(x, y_nodes, label = "f_b=$ang_freq", color = line_colors[i])
end

axislegend(position = :rt)
xlims!(0, 12)

CairoMakie.save(string(MA_DIR,"/rtsgmlc_inertia_vs_number_node_failures_show_title=$show_title.pdf"),fig)
CairoMakie.save(string(MA_DIR,"/rtsgmlc_inertia_vs_number_node_failures_show_title=$show_title.png"),fig)
fig
