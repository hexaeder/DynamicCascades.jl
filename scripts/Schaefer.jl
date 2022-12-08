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

# using GLMakie
# GLMakie.activate!()

# define pik colors
using Colors
Colors.color_names["pikorange"] = (227, 114, 34)
Colors.color_names["pikgray"] = (142, 144, 143)
Colors.color_names["pikcyan"] = (0, 159, 218)
Colors.color_names["pikgreen"] = (105, 146, 58)
Colors.color_names["pikblue"] = (124, 174, 175)

####
#### using two steadystates
####
network = import_system(:schaefer2018; γ = 0.8u"s")
sol1 = simulate(network;
    initial_fail = [5],
    failtime = 1.0,
    trip_lines = :none,
    tspan = (0.0, 40.0),
    solverargs = (; dtmax = 0.01));

sol2 = simulate(network;
    initial_fail = [5],
    failtime = 1.0,
    trip_lines = :dynamic,
    tspan = (0.0, 40.0),
    solverargs = (; dtmax = 0.01));
# inspect_solution(sol1)
t1 = 0.0
t2 = sol1.sol.t[end]

set_theme!(theme_minimal(), fontsize=20, textsize=25)
fig = Figure(resolution=(1500, 600))
fig[1,1] = leftpane = GridLayout()
node_size = [60, 70, 60, 60, 70]
# node_color = [colorant"pikblue", colorant"pikorange", colorant"pikblue", colorant"pikblue", colorant"pikorange"]
node_color = :black
edgecolortype = :relrating
gpargs = gparguments(sol1, t1; colortype = :relrating)
leftpane[1,1] = ax = Axis(fig)
p = graphplot!(ax, network; gpargs..., node_size, node_color,
               nlabels=repr.(1:5), nlabels_align=(:center, :center), nlabels_color=:white,
               nlabels_textsize=25)
ax.aspect = DataAspect()
p.edge_width[] = 2 .* p.edge_width[]
hidedecorations!(ax);
hidespines!(ax);
xlims!(ax, -1.05, 1.15)
ylims!(ax, -0.1, 1.5)
leftpane[2, 1] = Colorbar(fig, height = 25, vertical = false,
    colormap = DynamicCascades.edge_colorsheme(edgecolortype), label = "initial load relative to rating")
pos = GraphMakie.interpolate(get_edge_plot(p).paths[][5], 0.5)
scatter!(ax, pos, marker = '𐄂', color = :black, markersize = 100)

function plot_flows(ax, sol)
    for (idx,e) in enumerate(edges(sol.network))
        (ts, S) = seriesforidx(sol.load_S, idx)
        S .= S./0.978
        lines!(ax, ts, S, linewidth = 5, label=repr(e))
        xlims!(ax, (0, 15))
        hlines!(ax, 1.0, color=:black)
    end
end
fig[1,2] = rightpane = GridLayout()
rightpane[1,1] = ax1 = Axis(fig, ylabel="relative load")
rightpane[2,1] = ax2 = Axis(fig, ylabel="relative load", xlabel="time t (s)")
plot_flows(ax1, sol1)
plot_flows(ax2, sol2)
axislegend(ax2, position=:rc)
ylims!(ax1, 0, 1.1)
ylims!(ax2, 0, 1.1)

save(joinpath(PLOT_DIR, "schaefer_network.pdf"), fig)
