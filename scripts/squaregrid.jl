using DynamicCascades
using Graphs
using MetaGraphs
using Statistics
using GraphMakie
using Colors
using ColorSchemes
using DynamicCascades: PLOT_DIR
using Unitful
using NetworkDynamics

# using CairoMakie

using GLMakie
GLMakie.activate!()
set_theme!(theme_minimal(), fontsize = 20)

# 20 => 370
# 30 => 855

size = 30
x = 30
y = 30
network = import_system(:square; x, y, M=0.5u"s^2", D=0.1u"s", K=-10u"pu", seed=2)

sol = simulate(network;
               initial_fail = [855],
               failtime = 0.1,
               trip_lines = :none,
               tspan = (0.0, 50.0),
               solverargs = (; dtmax = 0.01));
inspect_solution(sol)

#0.4 first max

function graphplot_axis(fig, t::Observable, Δω; title)
    ax = Axis(fig)
    edge_color = [:white for i in 1:ne(network)]
    edge_color[sol.initial_fail] .= :black
    gpargs = gparguments(sol, t;
                         Δω,
                         node_colorscheme = ColorSchemes.diverging_bwr_40_95_c42_n256)
    graphplot!(ax, network; gpargs...,
                           edge_color, node_size=12)
    ax.aspect = DataAspect()
    ax.title = title
    hidedecorations!(ax); hidespines!(ax)
    return ax
end

# using CairoMakie
# CairoMakie.activate!()
# GLMakie.activate!()
fig = Figure(resolution=(1500, 1800))
# t1 = Observable(0.57)
# t2 = Observable(1.0)
# t3 = Observable(1.5)
# t4 = Observable(3.0)
t1 = Observable(0.5)
t2 = Observable(1.5)
t3 = Observable(3.0)
t4 = Observable(50.0)
Δω = Observable(0.12)
fig[2,1] = nwax = graphplot_axis(fig, t1, Δω; title="Frequencies at t₁ = $(t1[])")
fig[2,2] = graphplot_axis(fig, t2, Δω; title="Frequencies at t₂ = $(t2[])")
fig[2,3] = graphplot_axis(fig, t3, Δω; title="Frequencies at t₃ = $(t3[])")
fig[2,4] = graphplot_axis(fig, t4, Δω; title="Frequencies at t₄ = $(t4[])")
fig[1,:] = Colorbar(fig, get_node_plot(nwax.scene.plots[2]), height=25,width=1300, vertical=false, label="Node frequency in rad/s",)
# nwax.height[] = 500

fig[3,:] = ax2 = Axis(fig,
                      title ="Frequency of selected nodes right of incidence",
                      xlabel="time t in s", ylabel="frequency ω in rad/s",
                      xticks = @lift(([0,4,5, $t1, $t2, $t3, $t4],
                                      ["0", "4","5", "t₁ = $(t1[])","t₂ = $(t2[])", "t₃ = $(t3[])", "t₄ = $(t4[])"])))
ax2.xticksvisible[]=true
ax2.yticksvisible[]=true

# ax2.height[] = 500
nd, = nd_model(network);
ωidx = idx_containing(nd, "ω");
selnodes = [436, 438, 440]
labels = ["node 1 to the right of incidence",
          "node 3 to the right of incidence",
          "node 5 to the right of incidence",]
for (n, l) in zip(selnodes, labels)
    idx = ωidx[n]
    series = seriesforidx(sol.sol, idx)
    lines!(ax2, series...; linewidth=5, label=l)
end
axislegend(ax2, position=:rt)
xlims!(ax2, 0, 5)
# vlines!(ax2, t1, color=:black, linewidth=2)
# vlines!(ax2, t2, color=:black, linewidth=2)
# vlines!(ax2, t3, color=:black, linewidth=2)

function graphplot_axis_edges(fig, t::Observable, scaling; title)
    ax = Axis(fig)
    gpargs = gparguments(sol, t;
                         ecolortype = :abssteady,
                         activeP=true,
                         ecolorscaling=scaling,
                         ecolorscheme = Observable(ColorScheme(ColorSchemes.diverging_bwr_40_95_c42_n256[129:end]))
                         )
    p = graphplot!(ax, network; gpargs...,
                   node_size=0,
                   edge_width=5)

    ax.aspect = DataAspect()
    ax.title = title
    hidedecorations!(ax); hidespines!(ax)
    return ax
end
scaling = Observable(0.22)
fig[4,1] = nwax2 = graphplot_axis_edges(fig, t1, scaling; title="Load difference at t₁ = $(t1[])")
fig[4,2] = graphplot_axis_edges(fig, t2, scaling; title="Load difference at t₂ = $(t2[])")
fig[4,3] = graphplot_axis_edges(fig, t3, scaling; title="Load difference at t₃ = $(t3[])")
fig[4,4] = graphplot_axis_edges(fig, t4, scaling; title="Load difference at t₄ = $(t4[])")
# fig[5,:] =
#     Colorbar(fig,
#              ColorScheme(ColorSchemes.diverging_bwr_40_95_c42_n256[128:end]),
#              height=25,width=1300, vertical=false, label="Node frequency in rad/s",)
# ep = get_edge_plot(nwax2.scene.plots[2])

save(joinpath(PLOT_DIR, "failure_propagation.pdf"), fig)
GLMakie.activate!()
CairoMakie.activate!()
