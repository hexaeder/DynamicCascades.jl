using Revise
using DynamicCascades
using Graphs
using MetaGraphs
using Statistics
using GraphMakie
using Colors
using DynamicCascades: PLOT_DIR
using NetworkDynamics
using Unitful
using Unitful: @u_str
using DataFrames
using CSV
using LinearAlgebra
using CairoMakie

init = 27
damping = 0.1u"s"
network = import_system(:rtsgmlc; damping, tconst = 0.01u"s")
sol = simulate(network;
    initial_fail = Int[init],
    trip_lines = :dynamic,
    trip_nodes = :dynamic,
    solverargs = (;dtmax=0.01), verbose = true);

times = [0.0, sol.failures_nodes.t..., 6.0]

function plotnetwork(fig, sol, t; line=nothing, offset=nothing, text=nothing, apos=.5, ecolortype=:abssteady,
                     ecolorscaling=Observable(1.0), offlinecolor=colorant"lightgray")
    ax = Axis(fig)
    hidedecorations!(ax), hidespines!(ax)
    xlims!(ax, -10, 7)
    ylims!(ax, -5, 7)
    gpargs = gparguments(sol, t isa Observable ? t : Observable(t);
                         # ecolortype=Observable(ecolortype),
                         Δω=Observable(3.0),
                         offlinewidth=5,
                         offlinecolor,
                         ecolorscaling,
                         node_size=15)
    p = graphplot!(ax, network; gpargs...)

    if line !==nothing
        pos = GraphMakie.interpolate(get_edge_plot(p).paths[][line], apos)
        offset = Point2(offset)
        gapA = 0.0
        gapB = 0.6
        posA = pos + gapA*normalize(offset)
        posT = pos + offset
        posB = posT - gapB*normalize(offset)
        lines!(ax, [posA, posB], color=:black, linewidth=2)
        text!(ax, text; position=posT, align=(:center, :center))
    end
    return ax
end

CairoMakie.activate!()
set_theme!(theme_minimal(), fontsize = 20)
fig = Figure(resolution=(1500,1500))

fig[2,1] = nwax = plotnetwork(fig, sol, times[1]; line=11, offset=(-4,-.5), text="(3)")
fig[2,1] = Label(fig, "(a) t = $(round(times[1],digits=2)) s", tellwidth=false, tellheight=false, halign=:left, valign=:top)

fig[2,2] = plotnetwork(fig, sol, times[2]; line=27, offset=(0.75,2.5), text="(1)")
fig[2,2] = Label(fig, "(b) t = $(round(times[2],digits=2)) s", tellwidth=false, tellheight=false, halign=:left, valign=:top)

fig[2,3] = plotnetwork(fig, sol, times[3]; line=29, offset=(-2,2.5), text="(2)",apos=.4)
fig[2,3] = Label(fig, "(c) t = $(round(times[3],digits=2)) s", tellwidth=false, tellheight=false, halign=:left, valign=:top)

fig[3,1] = plotnetwork(fig, sol, times[4]; line=11, offset=(-4,-.5), text="(3)")
fig[3,1] = Label(fig, "(d) t = $(round(times[4],digits=2)) s", tellwidth=false, tellheight=false, halign=:left, valign=:top)

fig[3,2] = plotnetwork(fig, sol, times[5]; line=12, offset=(-4,0), text="(4)", apos=.25)
fig[3,2] = Label(fig, "(e) t = $(round(times[5],digits=2)) s", tellwidth=false, tellheight=false, halign=:left, valign=:top)

fig[3,3] = plotnetwork(fig, sol, times[6]; line=37, offset=(1.75,1.75), text="(5)")
fig[3,3] = Label(fig, "(f) t = $(round(times[6],digits=2)) s", tellwidth=false, tellheight=false, halign=:left, valign=:top)

fig[4,1] = plotnetwork(fig, sol, times[end-1])
fig[4,1] = Label(fig, "(g) t → ∞ (steady state)", tellwidth=false, tellheight=false, halign=:left, valign=:top)


fig[0,:] = Colorbar(fig, height=25, vertical=false,
                     colormap=Makie.ColorScheme([colorant"yellow", colorant"red"]), label="line load relative to rating")

fig[1,:] = Colorbar(fig, height=25, vertical=false,
                     colormap=ColorSchemes.diverging_bkr_55_10_c35_n256, label="node frequency deviation in hz")

save(joinpath(PLOT_DIR, "rts_cascade_nodes3.pdf"), fig)
