using Revise
using DynamicCascades
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

# only line failures
directory = "/home/vollmich/.julia/dev/MA_data/results_NB/20230722_145359.186inertia_vs_line_failures"
df_all_failures = DataFrame(CSV.File(string(directory,"/all_failures.csv")))

# plot data
fig = Figure(fontsize = 20)
fig[1, 1] = ax = Axis(fig;
    # title = L"Decreasing $G_{av}$ for one sample grid",
    # titlesize = 30,
    xlabel = "inertia M [s^2]",
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
    xlabel = "inertia M [s^2]",
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
    xlabel = "inertia M [s^2]",
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
