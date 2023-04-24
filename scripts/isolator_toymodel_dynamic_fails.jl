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

# for interactive plots
# using GLMakie
# GLMakie.activate!()

####
#### using two steadystates
####
γ = 0.8u"s"
M = 1u"s^2"
network = import_system(:isolator_toymodel; γ, M)
init_fail = 5
sol1 = simulate(network;
    initial_fail = [init_fail],
    failtime = 1.0,
    trip_lines = :none,
    trip_nodes = :none,
    tspan = (0.0, 15.0),
    solverargs = (; dtmax = 0.01));

sol2 = simulate(network;
    initial_fail = [init_fail],
    failtime = 1.0,
    trip_lines = :dynamic,
    trip_nodes = :none,
    tspan = (0.0, 15.0),
    solverargs = (; dtmax = 0.01));
# display(inspect_solution(sol1))
t1 = 0.0

set_theme!(theme_minimal(), fontsize=20, textsize=25)
fig = Figure(resolution=(1500, 600))
fig[1,1] = leftpane = GridLayout()
node_size = [40, 50, 40, 40, 50, 40, 50, 40, 40, 50,]
# node_color = [colorant"pikblue", colorant"pikorange", colorant"pikblue", colorant"pikblue", colorant"pikorange"]
node_color = :black
edgecolortype = :relrating
# initial load
leftpane[1,1] = ax = Axis(fig)
# gpargs = gparguments(sol1, t1; colortype = :relrating)
gpargs = gparguments(sol2, t1; colortype = :relrating, ecolorscaling = Observable(1.0), offlinecolor=colorant"lightgray")
p = graphplot!(ax, network; gpargs..., node_size, node_color,
               nlabels=repr.(1:10), nlabels_align=(:center, :center), nlabels_color=:white,
               nlabels_textsize=18)
ax.aspect = DataAspect()
p.edge_width[] = 1 .* p.edge_width[]
hidedecorations!(ax);
hidespines!(ax);
xlims!(ax, -1.05, 4.05)
ylims!(ax, -0.1, 1.5)
# final load
leftpane[2,1] = ax_final = Axis(fig)
# gpargs_final = gparguments(sol2, t2; colortype = :relrating)
gpargs_final = gparguments(sol2, sol2.sol.t[end]; colortype = :relrating, ecolorscaling = Observable(1.0), offlinecolor=colorant"lightgray")
p = graphplot!(ax_final, network; gpargs_final..., node_size, node_color,
               nlabels=repr.(1:10), nlabels_align=(:center, :center), nlabels_color=:white,
               nlabels_textsize=18)
ax_final.aspect = DataAspect()
p.edge_width[] = 1 .* p.edge_width[]
hidedecorations!(ax_final);
hidespines!(ax_final);
xlims!(ax_final, -1.05, 4.05)
ylims!(ax_final, -0.1, 1.5)

leftpane[3, 1] = Colorbar(fig, height = 25, vertical = false,
    colormap = DynamicCascades.edge_colorsheme(edgecolortype), label = "DYNAMIC failures: initial load (top), final load (bottom) relative to rating")
pos = GraphMakie.interpolate(get_edge_plot(p).paths[][init_fail], 0.5)
scatter!(ax, pos, marker = '𐄂', color = :black, markersize = 100)
scatter!(ax_final, pos, marker = '𐄂', color = :black, markersize = 100)

function plot_flows(ax, sol)
    for (idx,e) in enumerate(edges(sol.network))
        (ts, S) = seriesforidx(sol.load_S, idx)
        # S .= S./0.978
        S .= S./0.872 # 0.872 = K * tolerance (rating)
        lines!(ax, ts, S, linewidth = 2, label=repr(e))
        # xlims!(ax, (0, 15))
        xlims!(ax, (0, 15))
        hlines!(ax, 1.0, color=:black)
    end
end
fig[1,2] = rightpane = GridLayout()
rightpane[1,1] = ax1 = Axis(fig, ylabel="relative load")
rightpane[2,1] = ax2 = Axis(fig, ylabel="relative load", xlabel="time t (s)")
plot_flows(ax1, sol1)
plot_flows(ax2, sol2)
axislegend(ax1, position=:rt)
ylims!(ax1, 0, 1.2)
ylims!(ax2, 0, 1.1)
save(joinpath(PLOT_DIR, "new2_isolator_toymodel_dynamic_fails_D=$γ M=$M.pdf"), fig)
