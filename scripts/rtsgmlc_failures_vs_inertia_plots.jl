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
    xlabel = "inertia M [s^2]",
    # xlabelsize = 30,
    ylabel = "normalized average of line failures",
    # ylabelsize = 30
)

x = df_inertia_vs_failures.scale_inertia_values
y = df_inertia_vs_failures.rel_failures

# scatter!(x, y, color = :blue,label = "Test")
scatter!(x, y, color = :blue)
# axislegend()

CairoMakie.save(string(MA_DIR,"/inertia_vs_number_line_failures.pdf"),fig)

############################### node failures only #############################






directories = ["/results_plotting/node_failures_only/20230729_034424.446inertia_vs_line+node_failures_f_bound=0.24"]

ang_freq_bounds = [1.5, 2.0, 2.5, 3.0]

for (dir, val) in zip(directories, ang_freq_bounds)

directory = string(MA_DATA, dir)


df_all_failures_nodes = DataFrame(CSV.File(string(directory,"/all_failures_nodes.csv")))

# calculate relative number of failures
rel_number_failures = Float64[]
network = import_system(:rtsgmlc; damping=0.1u"s", tconst = 0.01u"s")
scale_inertia_values = names(df_all_failures)
for scale_inertia in scale_inertia_values
    number_failures_nodes = df_all_failures_nodes[!, scale_inertia]
    # NOTE TODO `mean(number_failures_nodes)` is not correct as long as line 11 is not included
    push!(rel_number_failures, mean(number_failures_nodes)/length(gen_node_idxs))
end
df_inertia_vs_failures_nodes = DataFrame("scale_inertia_values" => scale_inertia_values, "rel_failures" => rel_number_failures)
CSV.write(string(directory,"/inertia_vs_failures_nodes.csv"), df_inertia_vs_failures_nodes)

# load data
df_inertia_vs_failures_nodes = DataFrame(CSV.File(string(directory,"/inertia_vs_failures_nodes.csv")))

# plot data
fig = Figure(fontsize = 30)
Axis(fig[1, 1],
    # title = L"Decreasing $G_{av}$ for one sample grid",
    # titlesize = 30,
    xlabel = "inertia M [s^2]",
    # xlabelsize = 30,
    ylabel = "normalized average of node failures",
    # ylabelsize = 30
)

x = df_inertia_vs_failures.scale_inertia_values
y_nodes = df_inertia_vs_failures_nodes.rel_failures

scatter!(x, y_nodes, label = "nodes")
axislegend(position = :rb)


CairoMakie.save(string(MA_DIR,"/inertia_vs_number_line+node_failures.pdf"),fig)



######################### plot singe line failures #############################
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
