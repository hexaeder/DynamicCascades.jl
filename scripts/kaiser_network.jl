using Serialization
using DynamicCascades
using GraphMakie
using LightGraphs
using MetaGraphs
using ProgressMeter
using DataFrames
using Statistics
using GLMakie

using CairoMakie
CairoMakie.activate!()
DIR = "/Users/hw/MA/Forschungsbeleg/figures/"
set_theme!(theme_minimal(), fontsize=20)
set_theme!(resolution=(800, 600))

network1 = import_system(:kaiser2020; gen_γ=0.8, load_τ=0.1, seed=4)
network2 = deepcopy(network1)
rem_edge!(network2, 2, 22)
rem_edge!(network2, 1, 21)
set_prop!(network2, 1, 22, :X, 0.005)
set_prop!(network2, 2, 21, :X, 0.005)

failedge = (2, 17)
# failedge = (22, 32)
# failedge = (2, 8)
failidx1 = edgeidx(network1, failedge...)
sol1 = simulate(network1;
               initial_fail=[failidx1],
               trip_lines=false,
               tspan=(0., 1000.),
               solverargs=(;));

failidx2 = edgeidx(network2, failedge...)
sol2 = simulate(network2;
               initial_fail=[failidx2],
               trip_lines=false,
               tspan=(0., 1000.),
               solverargs=(;));

colortype=Observable(:abssteady)
ecolorscaling=Observable(0.15)
fig, ax, p = graphplot(sol1, sol1.sol.t[end]; colortype, ecolorscaling)
ax.aspect=DataAspect()
hidedecorations!(ax); hidespines!(ax)
save(joinpath(DIR,"kaiser_static_rerout1.pdf"), fig)

fig, ax, p = graphplot(sol2, sol2.sol.t[end]; colortype, ecolorscaling)
ax.aspect=DataAspect()
hidedecorations!(ax); hidespines!(ax)
save(joinpath(DIR,"kaiser_static_rerout2.pdf"), fig)

fig, ax, p = graphplot(sol1, sol1.sol.t[begin]; colortype)
ax.aspect=DataAspect()
hidedecorations!(ax); hidespines!(ax)
save(joinpath(DIR,"kaiser_gray1.pdf"), fig)

fig, ax, p = graphplot(sol2, sol2.sol.t[begin]; colortype)
ax.aspect=DataAspect()
hidedecorations!(ax); hidespines!(ax)
save(joinpath(DIR,"kaiser_gray2.pdf"), fig)
