using Revise
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

############################### line failures only #############################
directory = string(MA_DATA, "/results_plotting/line_failures_only/20230727_233120.943inertia_vs_line_failures")
df_all_failures = DataFrame(CSV.File(string(directory,"/all_failures.csv")))

# calculate relative number of failures
rel_number_failures = Float64[]
network = import_system(:rtsgmlc; damping=0.1u"s", tconst = 0.01u"s")
scale_inertia_values = names(df_all_failures)
for scale_inertia in scale_inertia_values
    number_failures = df_all_failures[!, scale_inertia]
    push!(rel_number_failures, mean(number_failures)/(ne(network)-1))
end

df_inertia_vs_failures = DataFrame("scale_inertia_values" => scale_inertia_values, "rel_failures" => rel_number_failures)
CSV.write(string(directory,"/inertia_vs_failures.csv"), df_inertia_vs_failures)

# load data
df_inertia_vs_failures = DataFrame(CSV.File(string(directory,"/inertia_vs_failures.csv")))

# plot data
fig = Figure(fontsize = 30)
Axis(fig[1, 1],
    # title = L"Decreasing $G_{av}$ for one sample grid",
    # titlesize = 30,
    xlabel = "scaling factor of inertia",
    # xlabelsize = 30,
    ylabel = "normalized average of line failures",
    # ylabelsize = 30
)

x = df_inertia_vs_failures.scale_inertia_values
y = df_inertia_vs_failures.rel_failures

# scatter!(x, y, color = :blue, label = "Test")
scatterlines!(x, y, color = :blue, linewidth=2)
# axislegend()

CairoMakie.save(string(MA_DIR,"/rtsgmlc_inertia_vs_number_line_failures.pdf"),fig)

############################### node failures only #############################
directories = ["/results_plotting/node_failures_only/20230729_034424.446inertia_vs_line+node_failures_f_bound=0.24",
            "/results_plotting/node_failures_only/20230729_034425.247inertia_vs_line+node_failures_f_bound=0.32",
            "/results_plotting/node_failures_only/20230729_034431.966inertia_vs_line+node_failures_f_bound=0.4",
            "/results_plotting/node_failures_only/20230729_034435.458inertia_vs_line+node_failures_f_bound=0.48",
            "/results_plotting/node_failures_only/20230730_120321.717inertia_vs_line+node_failures_f_bound=0.64",
            "/results_plotting/node_failures_only/20230730_121128.261inertia_vs_line+node_failures_f_bound=0.8",
            "/results_plotting/node_failures_only/20230730_121128.342inertia_vs_line+node_failures_f_bound=0.95",
            "/results_plotting/node_failures_only/20230730_121128.333inertia_vs_line+node_failures_f_bound=1.27",
            "/results_plotting/node_failures_only/20230730_121128.315inertia_vs_line+node_failures_f_bound=1.59"
            # "/results_plotting/node_failures_only/20230730_121159.099inertia_vs_line+node_failures_f_bound=1.91"
            ]

ang_freq_bounds = [0.24, 0.32, 0.40, 0.48, 0.64, 0.80, 0.95, 1.27, 1.59, 1.91]

norm_values = (ang_freq_bounds .- minimum(ang_freq_bounds)) ./ (maximum(ang_freq_bounds) - minimum(ang_freq_bounds))

# Function to generate distinct colors based on a color map and normalization
function distinct_colors(color_map, values)
    norm_values = (values .- minimum(values)) ./ (maximum(values) - minimum(values))

    colors = [cgrad(color_map, 101; categorical = true, rev=true)[Int(ceil(i*100)+1)] for i in norm_values]
    return colors
end


# color_map = ColorSchemes.plasma
color_map = ColorSchemes.cividis
# color_map = :blues

# Generate distinct colors based on the ang_freq_bounds
line_colors = distinct_colors(color_map, ang_freq_bounds)

# plot data
fig = Figure(fontsize = 28)
Axis(fig[1, 1],
    title = "Variation of angular frequency bounds",
    # titlesize = 30,
    xlabel = "scaling factor of inertia",
    # xlabelsize = 30,
    ylabel = "normalized average of node failures",
    # ylabelsize = 30
)

network = import_system(:rtsgmlc; damping = 0.1u"s", tconst = 0.01u"s")
nd, = nd_model(network)
ω_state_idxs = idx_containing(nd, "ω")
gen_node_idxs = map(s -> parse(Int, String(s)[4:end]), nd.syms[ω_state_idxs])
lines!([NaN], [NaN]; label="Bounds", color=:white, linewidth=3)

for (i, dir, ang_freq) in zip(1:length(line_colors),directories, ang_freq_bounds)
    directory = string(MA_DATA, dir)
    df_all_failures_nodes = DataFrame(CSV.File(string(directory,"/all_failures_nodes.csv")))

    # calculate relative number of failures
    rel_number_failures = Float64[]
    scale_inertia_values = names(df_all_failures_nodes)
    for scale_inertia in scale_inertia_values
        number_failures_nodes = df_all_failures_nodes[!, scale_inertia]
        # NOTE TODO `mean(number_failures_nodes)` is not correct as long as line 11 is not included
        push!(rel_number_failures, mean(number_failures_nodes)/length(gen_node_idxs))
    end
    df_inertia_vs_failures_nodes = DataFrame("scale_inertia_values" => scale_inertia_values, "rel_failures" => rel_number_failures)
    CSV.write(string(directory,"/inertia_vs_failures_nodes.csv"), df_inertia_vs_failures_nodes)

    # load data
    df_inertia_vs_failures_nodes = DataFrame(CSV.File(string(directory,"/inertia_vs_failures_nodes.csv")))

    x = df_inertia_vs_failures_nodes.scale_inertia_values
    y_nodes = df_inertia_vs_failures_nodes.rel_failures

    scatterlines!(x, y_nodes, label = "$ang_freq", color = line_colors[i])
end

axislegend(position = :rt)
# ylims!(-0.0005, 0.002)
xlims!(0, 12)
fig

CairoMakie.save(string(MA_DIR,"/rtsgmlc_inertia_vs_number_node_failures.pdf"),fig)

############################### node + line failures ###########################

# all bounds
directories = [
            "/results_plotting/line+node_failures/20230728_013046.031inertia_vs_line+node_failures_f_bound=0.16",
            "/results_plotting/line+node_failures/20230728_013105.59inertia_vs_line+node_failures_f_bound=0.32",
            "/results_plotting/line+node_failures/20230728_013105.412inertia_vs_line+node_failures_f_bound=0.48",
            "/results_plotting/line+node_failures/20230730_082824.996inertia_vs_line+node_failures_f_bound=0.64",
            "/results_plotting/line+node_failures/20230730_083925.633inertia_vs_line+node_failures_f_bound=0.8",
            "/results_plotting/line+node_failures/20230730_084103.401inertia_vs_line+node_failures_f_bound=0.95",
            "/results_plotting/line+node_failures/20230730_083925.569inertia_vs_line+node_failures_f_bound=1.27",
            "/results_plotting/line+node_failures/20230730_083925.646inertia_vs_line+node_failures_f_bound=1.59",
            "/results_plotting/line+node_failures/20230730_083955.126inertia_vs_line+node_failures_f_bound=1.91"
            ]
ang_freq_bounds = [0.16, 0.32, 0.48, 0.64, 0.80, 0.95, 1.27, 1.59, 1.91]

# lower bounds
directories = [
            "/results_plotting/line+node_failures/20230728_013105.412inertia_vs_line+node_failures_f_bound=0.48"
            ]
ang_freq_bounds = [0.48]

# lower bounds
directories = [
            "/results_plotting/line+node_failures/20230728_013046.031inertia_vs_line+node_failures_f_bound=0.16",
            "/results_plotting/line+node_failures/20230728_013105.59inertia_vs_line+node_failures_f_bound=0.32",
            "/results_plotting/line+node_failures/20230728_013105.412inertia_vs_line+node_failures_f_bound=0.48"
            ]
ang_freq_bounds = [0.16, 0.32, 0.48]

# larger bounds
directories = [
            # "/results_plotting/line+node_failures/20230730_082824.996inertia_vs_line+node_failures_f_bound=0.64",
            # "/results_plotting/line+node_failures/20230730_083925.633inertia_vs_line+node_failures_f_bound=0.8",
            "/results_plotting/line+node_failures/20230730_084103.401inertia_vs_line+node_failures_f_bound=0.95",
            "/results_plotting/line+node_failures/20230730_083925.569inertia_vs_line+node_failures_f_bound=1.27",
            # "/results_plotting/line+node_failures/20230730_083925.646inertia_vs_line+node_failures_f_bound=1.59",
            "/results_plotting/line+node_failures/20230730_083955.126inertia_vs_line+node_failures_f_bound=1.91"
            ]
ang_freq_bounds = [0.64, 0.95, 1.27,  1.91]

# only one wide and one narrow bound
directories = [
            # "/results_plotting/line+node_failures/20230728_013105.59inertia_vs_line+node_failures_f_bound=0.32",
            "/results_plotting/line+node_failures/20230728_013105.412inertia_vs_line+node_failures_f_bound=0.48",
            "/results_plotting/line+node_failures/20230730_082824.996inertia_vs_line+node_failures_f_bound=0.64"
            # "/results_plotting/line+node_failures/20230730_083925.633inertia_vs_line+node_failures_f_bound=0.8"
            ]
ang_freq_bounds = [0.48, 0.64]


norm_values = (ang_freq_bounds .- minimum(ang_freq_bounds)) ./ (maximum(ang_freq_bounds) - minimum(ang_freq_bounds))

# color_map = ColorSchemes.plasma
# color_map = ColorSchemes.cividis
# color_map = :blues

# Generate distinct colors based on the ang_freq_bounds
# line_colors = distinct_colors(color_map, ang_freq_bounds)
line_colors = [Makie.wong_colors()[2], Makie.wong_colors()[3]]

# NOTE uncomment for 0.48
# line_colors =[:blue]

# plot data
fig = Figure(fontsize = 28)
Axis(fig[1, 1],
    title = "Summed line and node failures",
    # titlesize = 30,
    xlabel = "scaling factor of inertia",
    # xlabelsize = 30,
    ylabel = "average of failures",
    # ylabelsize = 30
)


network = import_system(:rtsgmlc; damping = 0.1u"s", tconst = 0.01u"s")
nd, = nd_model(network)
ω_state_idxs = idx_containing(nd, "ω")
gen_node_idxs = map(s -> parse(Int, String(s)[4:end]), nd.syms[ω_state_idxs])
# lines!([NaN], [NaN]; label="Bounds", color=:white, linewidth=3)

for (i, dir, ang_freq) in zip(1:length(line_colors),directories, ang_freq_bounds)
    directory = string(MA_DATA, dir)

    # calculate relative number of failures
    # nodes
    df_all_failures_nodes = DataFrame(CSV.File(string(directory,"/all_failures_nodes.csv")))
    rel_number_failures_nodes = Float64[] # NOTE this is actually the mean
    scale_inertia_values = names(df_all_failures_nodes)
    for scale_inertia in scale_inertia_values
        number_failures_nodes = df_all_failures_nodes[!, scale_inertia]
        # NOTE TODO `mean(number_failures_nodes)` is not correct as long as line 11 is not included
        # # normalized
        # push!(rel_number_failures_nodes, mean(number_failures_nodes)/length(gen_node_idxs))
        # not normalized
        push!(rel_number_failures_nodes, mean(number_failures_nodes))
    end
    df_inertia_vs_failures_nodes = DataFrame("scale_inertia_values" => scale_inertia_values, "rel_failures" => rel_number_failures_nodes)
    CSV.write(string(directory,"/inertia_vs_failures_nodes.csv"), df_inertia_vs_failures_nodes)

    #lines
    df_all_failures_lines = DataFrame(CSV.File(string(directory,"/all_failures.csv")))
    rel_number_failures_lines = Float64[] # NOTE this is actually the mean
    for scale_inertia in scale_inertia_values
        number_failures_lines = df_all_failures_lines[!, string(scale_inertia)]
        # # normalized
        # push!(rel_number_failures_lines, mean(number_failures_lines)/(ne(network)-1))
        # not normalized
        push!(rel_number_failures_lines, mean(number_failures_lines))
    end
    df_inertia_vs_failures_lines = DataFrame("scale_inertia_values" => scale_inertia_values, "rel_failures" => rel_number_failures_lines)
    CSV.write(string(directory,"/inertia_vs_failures.csv"), df_inertia_vs_failures_lines)

    # load data
    df_inertia_vs_failures_nodes = DataFrame(CSV.File(string(directory,"/inertia_vs_failures_nodes.csv")))
    df_inertia_vs_failures = DataFrame(CSV.File(string(directory,"/inertia_vs_failures.csv")))

    x = df_inertia_vs_failures_nodes.scale_inertia_values
    y_nodes = df_inertia_vs_failures_nodes.rel_failures
    y_lines = df_inertia_vs_failures_lines.rel_failures

    scatterlines!(x, (y_nodes + y_lines), label = "f_b=$ang_freq", color = line_colors[i], linewidth=4, markersize=10)
    # scatterlines!(x, y_lines, color = line_colors[i], linestyle=:dash )
    # scatterlines!(x, y_nodes, label = "$ang_freq", color = line_colors[i])
    # scatterlines!(x, y_lines, color = line_colors[i], linestyle=:dash )
end
# lines!([NaN], [NaN]; label="Lines dashed", color=:black, linewidth=3, linestyle=:dash)

# all bounds
# text!(9.0, 0.5, text = "node failures: solid lines ___ \n line failures: dashed lines -----", align = (:center, :center), textsize=25)
# axislegend(position = :rt)

# all bounds zoomed in
# text!(9.0, 0.08, text = "node failures: solid lines ___ \n line failures: dashed lines -----", align = (:center, :center), textsize=25)
axislegend(position = :rc)
# ylims!(-0.0005, 0.1)

# low bounds
# text!(16.0, 0.3, text = "node failures: solid lines ___ \n line failures: dashed lines -----", align = (:center, :center), textsize=25)
# axislegend(position = :rt)

# large bounds
# text!(9.0, 0.004, text = "node failures: solid lines ___ \n line failures: dashed lines -----", align = (:center, :center), textsize=25)
# axislegend(position = :rb)

# 0.48
# text!(15.0, 0.01, text = "node failures: solid lines ___ \n line failures: dashed lines -----", align = (:center, :center), textsize=25)
# axislegend(position = :rt)

# ylims!(-0.0005, 0.1)
# xlims!(0, 5)
fig


# CairoMakie.save(string(MA_DIR,"/rtsgmlc_inertia_vs_number_line+node_failures_low_bounds.pdf"),fig)
# CairoMakie.save(string(MA_DIR,"/rtsgmlc_inertia_vs_number_line+node_failures_large_bounds.pdf"),fig)
# CairoMakie.save(string(MA_DIR,"/rtsgmlc_inertia_vs_number_line+node_failures_0.48.pdf"),fig)
# CairoMakie.save(string(MA_DIR,"/rtsgmlc_inertia_vs_number_line+node_failures_all_bounds.pdf"),fig)
# CairoMakie.save(string(MA_DIR,"/rtsgmlc_inertia_vs_number_line+node_failures_all_bounds_zoomed_in.pdf"),fig)
CairoMakie.save(joinpath(MA_DIR, "rtsgmlc_uebergang_inertia_line+node_failures.pdf"),fig)
CairoMakie.save(joinpath(MA_DIR, "rtsgmlc_uebergang_inertia_line+node_failures.png"),fig)

######################### plot singe line failures #############################
directory = "/home/vollmich/.julia/dev/MA_data/results_NB/20230722_145359.186inertia_vs_line_failures"
df_all_failures = DataFrame(CSV.File(string(directory,"/all_failures.csv")))

# plot data
fig = Figure(fontsize = 20)
fig[1, 1] = ax = Axis(fig;
    # title = L"Decreasing $G_{av}$ for one sample grid",
    # titlesize = 30,
    xlabel = "scaling factor of inertia",
    xlabelsize = 30,
    ylabel = "normalized number of failures",
    ylabelsize = 30
)

num_rows = size(df_all_failures, 1)
for i in 1:num_rows
    if sum(values(df_all_failures[i,:])) !== 0.0
        y_lines = 1/num_rows .* Vector(df_all_failures[i,:])
        lines!(ax, x, y_lines; label="line $i", fontsize=10)
    end
end

axislegend()
fig

save(joinpath(MA_DIR, "rtsgmlc_failures_vs_inertia_single_lines_only.pdf"), fig)
save(joinpath(MA_DIR, "rtsgmlc_failures_vs_inertia_single_lines_only.png"), fig)

# line and node failures
# line failures
directory = "/home/vollmich/.julia/dev/MA_data/results_NB/20230724_232058.139inertia_vs_line+node_failures"
df_all_failures = DataFrame(CSV.File(string(directory,"/all_failures.csv")))

# plot data
fig = Figure(fontsize = 20)
fig[1, 1] = ax = Axis(fig;
    # title = L"Decreasing $G_{av}$ for one sample grid",
    # titlesize = 30,
    xlabel = "scaling factor of inertia",
    xlabelsize = 30,
    ylabel = "normalized number of line failures",
    ylabelsize = 30
)

num_rows = size(df_all_failures, 1)
for i in 1:num_rows
    if sum(values(df_all_failures[i,:])) !== 0.0
        y_lines = 1/num_rows .* Vector(df_all_failures[i,:])
        lines!(ax, x, y_lines; label="line $i", fontsize=10)
    end
end

axislegend()
fig
save(joinpath(MA_DIR, "rtsgmlc_failures_vs_inertia_single_lines+nodes_lines.pdf"), fig)
save(joinpath(MA_DIR, "rtsgmlc_failures_vs_inertia_single_lines+nodes_lines.png"), fig)

# node failures
df_all_failures = DataFrame(CSV.File(string(directory,"/all_failures_nodes.csv")))

# plot data
fig = Figure(fontsize = 20)
fig[1, 1] = ax = Axis(fig;
    # title = L"Decreasing $G_{av}$ for one sample grid",
    # titlesize = 30,
    xlabel = "scaling factor of inertia",
    xlabelsize = 30,
    ylabel = "normalized number of node failures",
    ylabelsize = 30
)

damping = 0.1u"s"
network = import_system(:rtsgmlc; damping, tconst = 0.01u"s")
nd, = nd_model(network)
ω_state_idxs = idx_containing(nd, "ω")
gen_node_idxs = map(s -> parse(Int, String(s)[4:end]), nd.syms[ω_state_idxs])
num_gen_nodes = length(gen_node_idxs)
num_rows = size(df_all_failures, 1)
for i in 1:num_rows
    if sum(values(df_all_failures[i,:])) !== 0.0
        y_lines = 1/num_gen_nodes .* Vector(df_all_failures[i,:])
        lines!(ax, x, y_lines; label="line $i", fontsize=10)
    end
end

axislegend()
fig
save(joinpath(MA_DIR, "rtsgmlc_failures_vs_inertia_single_lines+nodes_nodes.pdf"), fig)
save(joinpath(MA_DIR, "rtsgmlc_failures_vs_inertia_single_lines+nodes_nodes.png"), fig)
