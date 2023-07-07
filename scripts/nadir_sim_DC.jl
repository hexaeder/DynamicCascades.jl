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

####
#### using two steadystates
####
network = import_system(:nadir_sim; M=5.0u"s^2", γ=0.1u"s", tconst=0.01u"s")
sol = simulate(network;
    initial_fail = [1],
    failtime = 1.0,
    trip_lines = :none,
    trip_nodes = :none,
    trip_load_nodes = :none,
    tspan = (0.0, 200.0),
    solverargs = (; dtmax = 0.01));

# plot solution
function plotnetwork(fig, sol, t)
    ax = Axis(fig)
    hidedecorations!(ax), hidespines!(ax)
    xlims!(ax, -10, 7)
    ylims!(ax, -5, 7)
    gpargs = gparguments(sol, t;
                         colortype=:abssteady,
                         Δω=Observable(2.5),
                         offlinewidth=3,
                         offlinecolor=colorant"lightgray",
                         ecolorscaling = Observable(1.0),
                         node_size=15,
                         show_labels=false)
    p = graphplot!(ax, network; gpargs...)
    return ax, p
end

# tobs = Observable(sol.sol.t[end-13])
# tobs = Observable(sol.sol.t[end])
tobs = Observable(0.0)
fig = Figure(resolution=(1800,1000))

# plot network
# fig[1,1] = Label(fig, @lift("t = "*repr(round($tobs,digits=2))*" s"), tellwidth=false)
fig[1,1], p = plotnetwork(fig, sol, tobs)

# # plot all power flows
# fig[2,1] = ax = Axis(fig; xlabel="time t in s", ylabel="apparent power flow in p.u.", title="power flow of all lines")
# for i in 1:length(sol.load_S.saveval[1])
#     t = sol.load_S.t
#     # t = sol.frequencies_load_nodes.t[1:20]
#     y = [sol.load_S.saveval[t][i] for t in 1:length(sol.load_S.t)]
#     # y = [sol.frequencies_load_nodes.saveval[t][i] for t in 1:20]
#     # lines!(ax, t, y; label="frequency on load node ($i)", linewidth=3)
#     scatter!(ax, t, y; label="power flow on line ($i)", markersize=5)
# end
# vlines!(ax, tobs; color=:black, linewidth=1)
# fig

# plot failed power flows
fig[2,1] = ax = Axis(fig; xlabel="time t in s", ylabel="apparent power flow in p.u.", title=@lift("t = "*repr(round($tobs,digits=2))*" s                                      flow transients of failing lines"))
for (i, l) in pairs(sol.failures.saveval)
    t, s = seriesforidx(sol.load_S, l)
    lines!(ax, t, s; label="flow on edge ($i)", linewidth=3)
    # scatter!(ax, t, s; label="flow on edge ($i)", markersize=5)
    scatter!(ax, (t[end], s[end]); marker=:star5, markersize=25)
end
vlines!(ax, tobs; color=:black, linewidth=1)

# plot frequencies of all gen nodes
fig[1,2] = ax = Axis(fig; xlabel="time t in s", ylabel="angular frequency ω", title="frequency transients of all generator nodes")
(nd, p, overload_cb) = nd_model(network)
state_idx = idx_containing(nd, "ω") # array: indices of ω-states
for i in state_idx
    t = sol.sol.t
    # t = sol.sol.t[1:300]
    y = [sol.sol.u[t][i] for t in 1:length(sol.sol.t)]
    # y = [sol.sol.u[t][i] for t in 1:300]
    lines!(ax, t, y; label="frequency on node ($i)", linewidth=3)
    # scatter!(ax, t, y; label="frequency on node ($i)", linewidth=3)
end
vlines!(ax, tobs; color=:black, linewidth=1)

# plot frequencies of failed gen nodes
fig[2,2] = ax = Axis(fig; xlabel="time t in s", ylabel="angular frequency ω", title="frequency transients of failing generator nodes")
(nd, p, overload_cb) = nd_model(network)
state_idx = idx_containing(nd, "ω") # array: indices of ω-states
node_idx = map(s -> parse(Int, String(s)[4:end]), nd.syms[state_idx]) # array: indices of vertices that are generators
for (i, l) in pairs([state_idx[findfirst(x -> x == i, node_idx)] for i in sol.failures_nodes.saveval])
    t, s = seriesforidx(sol.sol, l)
    lines!(ax, t, s; label="frequency on node ($i)", linewidth=3)
    scatter!(ax, (t[end], s[end]); marker=:star5, markersize=25)
end
vlines!(ax, tobs; color=:black, linewidth=1)

# plot frequencies of all load nodes
fig[1,3] = ax = Axis(fig; xlabel="time t in s", ylabel="angular frequency ω", title="frequency transients of all load nodes")
for i in 1:length(sol.frequencies_load_nodes.saveval[1])
    t = sol.frequencies_load_nodes.t
    # t = sol.frequencies_load_nodes.t[1:20]
    y = [sol.frequencies_load_nodes.saveval[t][i] for t in 1:length(sol.frequencies_load_nodes.t)]
    # y = [sol.frequencies_load_nodes.saveval[t][i] for t in 1:20]
    lines!(ax, t, y; label="frequency on load node ($i)", linewidth=3)
    # scatter!(ax, t, y; label="frequency on load node ($i)", markersize=5)
end
vlines!(ax, tobs; color=:black, linewidth=1)

# plot frequencies of failed load nodes
load_node_idxs = findall(x -> x==:load, get_prop(network, 1:nv(network), :type))
fig[2,3] = ax = Axis(fig; xlabel="time t in s", ylabel="angular frequency ω", title="frequency transients of failing load nodes")
for (i, l) in pairs([findfirst(x -> x == i, load_node_idxs) for i in sol.failures_load_nodes.saveval])
    t, s = seriesforidx(sol.frequencies_load_nodes, l)
    # scatter!(ax, t, s; label="frequency on load node ($i)", linewidth=3)
    lines!(ax, t, s; label="frequency on load node ($i)", linewidth=3)
    scatter!(ax, (t[end], s[end]); marker=:star5, markersize=25)
end
vlines!(ax, tobs; color=:black, linewidth=1)
fig




# display(inspect_solution(sol1))
#
# sol2 = simulate(network;
#     initial_fail = [3],
#     failtime = 1.0,
#     trip_lines = :dynamic,
#     trip_nodes = :none,
#     trip_load_nodes = :none,
#     tspan = (0.0, 40.0),
#     solverargs = (; dtmax = 0.01));
#
# t1 = 0.0
# t2 = sol1.sol.t[end]
#
#
#
# set_theme!(theme_minimal(), fontsize=20, textsize=25)
# fig = Figure(resolution=(1500, 600))
# fig[1,1] = leftpane = GridLayout()
# node_size = [90, 60]
# # node_color = [colorant"pikblue", colorant"pikorange", colorant"pikblue", colorant"pikblue", colorant"pikorange"]
# node_color = :black
# edgecolortype = :relrating
# gpargs = gparguments(sol1, t1; colortype = :relrating, ecolorscaling = Observable(1.0))
# leftpane[1,1] = ax = Axis(fig)
# p = graphplot!(ax, network; gpargs..., node_size, node_color,
#                nlabels=repr.(1:2), nlabels_align=(:center, :center), nlabels_color=:white,
#                nlabels_textsize=25)
# ax.aspect = DataAspect()
# p.edge_width[] = 2 .* p.edge_width[]
# hidedecorations!(ax);
# hidespines!(ax);
# xlims!(ax, -1.05, 1.15)
# ylims!(ax, -0.1, 1.5)
# leftpane[2, 1] = Colorbar(fig, height = 25, vertical = false,
#     colormap = DynamicCascades.edge_colorsheme(edgecolortype), label = "initial load relative to rating")
# # pos = GraphMakie.interpolate(get_edge_plot(p).paths[][4], 0.5)
# # scatter!(ax, pos, marker = '𐄂', color = :black, markersize = 100)
#
# function plot_flows(ax, sol)
#     for (idx,e) in enumerate(edges(sol.network))
#         (ts, S) = seriesforidx(sol.load_S, idx)
#         lines!(ax, ts, S, linewidth = 5, label=repr(e))
#         xlims!(ax, (0, 15))
#         hlines!(ax, 1.0, color=:black)
#     end
# end
# fig[1,2] = rightpane = GridLayout()
# rightpane[1,1] = ax1 = Axis(fig, ylabel="load")
# rightpane[2,1] = ax2 = Axis(fig, ylabel="load", xlabel="time t (s)")
# plot_flows(ax1, sol1)
# plot_flows(ax2, sol1)
# axislegend(ax2, position=:rc)
# # ylims!(ax1, 0, 1.1)
# # ylims!(ax2, 0, 1.1)
# xlims!(ax1, 0.9999, 1.0001)
# xlims!(ax2, 0.9999, 1.0001)
# fig
#
# save(joinpath(PLOT_DIR, "nadir_sim.pdf"), fig)
