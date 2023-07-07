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

using GLMakie
GLMakie.activate!()

####
#### using two steadystates
####
network = import_system(:toymodel_square; M=1.0u"s^2", γ=1.1u"s", tconst=0.01u"s")
sol1 = simulate(network;
    initial_fail = [4],
    failtime = 1.0,
    trip_lines = :none,
    trip_nodes = :none,
    trip_load_nodes = :none,
    tspan = (0.0, 200.0),
    solverargs = (; dtmax = 0.01));

display(inspect_solution(sol1))

sol2 = simulate(network;
    initial_fail = [3],
    failtime = 1.0,
    trip_lines = :dynamic,
    trip_nodes = :none,
    trip_load_nodes = :none,
    tspan = (0.0, 40.0),
    solverargs = (; dtmax = 0.01));

t1 = 0.0
t2 = sol1.sol.t[end]

######################## testing initial steady state ##########################
(nd, p, overload_cb) = nd_model(network)
nd.syms

x_static=steadystate(network)
δ1 = x_static[1]
ω1 = x_static[2]
δ2 = x_static[3]
ω1 = x_static[4]
δ3 = x_static[5]
ω1 = x_static[6]
δ4 = x_static[7]
ω1 = x_static[8]

ω_dot1 = +1 -sin(δ1-δ2) - sin(δ1-δ3)
ω_dot2 = -1 -sin(δ2-δ1) - sin(δ2-δ4)
ω_dot3 = -1 -sin(δ3-δ1) - sin(δ3-δ4)
ω_dot4 = +1 -sin(δ4-δ3) - sin(δ4-δ2)
################################################################################

set_theme!(theme_minimal(), fontsize=20, textsize=25)
fig = Figure(resolution=(1500, 600))
fig[1,1] = leftpane = GridLayout()
node_size = [90, 60, 60, 90]
# node_color = [colorant"pikblue", colorant"pikorange", colorant"pikblue", colorant"pikblue", colorant"pikorange"]
node_color = :black
edgecolortype = :relrating
gpargs = gparguments(sol1, t1; colortype = :relrating, ecolorscaling = Observable(1.0))
leftpane[1,1] = ax = Axis(fig)
p = graphplot!(ax, network; gpargs..., node_size, node_color,
               nlabels=repr.(1:4), nlabels_align=(:center, :center), nlabels_color=:white,
               nlabels_textsize=25)
ax.aspect = DataAspect()
p.edge_width[] = 2 .* p.edge_width[]
hidedecorations!(ax);
hidespines!(ax);
xlims!(ax, -1.05, 1.15)
ylims!(ax, -0.1, 1.5)
leftpane[2, 1] = Colorbar(fig, height = 25, vertical = false,
    colormap = DynamicCascades.edge_colorsheme(edgecolortype), label = "initial load relative to rating")
pos = GraphMakie.interpolate(get_edge_plot(p).paths[][4], 0.5)
scatter!(ax, pos, marker = '𐄂', color = :black, markersize = 100)

function plot_flows(ax, sol)
    for (idx,e) in enumerate(edges(sol.network))
        (ts, S) = seriesforidx(sol.load_S, idx)
        lines!(ax, ts, S, linewidth = 5, label=repr(e))
        xlims!(ax, (0, 15))
        hlines!(ax, 1.0, color=:black)
    end
end
fig[1,2] = rightpane = GridLayout()
rightpane[1,1] = ax1 = Axis(fig, ylabel="load")
rightpane[2,1] = ax2 = Axis(fig, ylabel="load", xlabel="time t (s)")
plot_flows(ax1, sol1)
plot_flows(ax2, sol1)
axislegend(ax2, position=:rc)
ylims!(ax1, 0, 1.1)
ylims!(ax2, 0, 1.1)
fig

save(joinpath(PLOT_DIR, "toymodel.pdf"), fig)
