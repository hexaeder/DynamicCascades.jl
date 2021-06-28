#!/usr/bin/env julia

using OrdinaryDiffEq
using DiffEqCallbacks
using SteadyStateDiffEq
using DataFrames

using DynamicCascades
using NetworkLayout
using LightGraphs
using MetaGraphs
using GLMakie
using GraphMakie

DIR = "/Users/hw/MA/Forschungsbeleg/figures/"
set_theme!(theme_minimal(), fontsize=20)

network = import_system(:rts96prepared; gen_γ=0.8, slack_γ=0.8, load_τ=0.1, losses=false)

sol = simulate(network;
               tspan=(0.0,500.0),
               initial_fail=[27],
               trip_lines=false);

x_static = sol.sol[end]
(nd, p) = nd_model(network)
p[2][27] = 0.0
S = calculate_apparent_power(x_static, p, 0.0, nd.f, network)
margin = describe_edges(network).rating .- S
@assert minimum(margin) > 0

fig = Figure(resolution=(1000,1000))
fig[1,1] = ax = Axis(fig)
hidedecorations!(ax), hidespines!(ax)
xlims!(-10, 7)
ylims!(-5, 7)

t = Node(1.0)
ecolorscaling = Node(0.5)
gpargs = gparguments(sol, t;
                     colortype=:abssteady,
                     offlinewidth=0.0,
                     ecolorscaling)
p = graphplot!(ax, network; gpargs...)
timelabel = @lift "time = " * repr(round($t, digits=2)) * " s"
fig[1,1] = Label(fig, timelabel, textsize=30,
                 halign=:left, valign=:top, tellwidth=false, tellheight=false)

function getline(idx)
    e = collect(edges(network))[idx]
    p[:node_positions][][[e.src, e.dst]]
end

largs = (;color=:black, linestyle=:dot, linewidth=3.0)
l27 = lines!(ax, getline(27); largs...)

file = "rts_wo_fail"
tspan = range(1.0, 30, length=10*30)
save(joinpath(DIR, "../videos", file*".png"), fig)
record(fig, joinpath(DIR, "../videos", file*".mp4"), tspan; framerate=30) do time
    t[] = time
end
save(joinpath(DIR, "../videos", file*"_end"*".png"), fig)

####
#### now for the actual simulation
####
sol = simulate(network;
               tspan=(0.0,500.0),
               initial_fail=[27],
               trip_lines=true);
# inspect_solution(sol)

times = [1.0, sol.failures.t[1:3]..., 50]
fps = 30
timeranges = []
for i in 1:length(times)-1
    tstart = times[i]
    tend = times[i+1]
    r = range(tstart, tend, length=round(Int, (tend-tstart)*fps))
    push!(timeranges, r)
end

fig = Figure(resolution=(1000,1000))
fig[1,1] = ax = Axis(fig)
hidedecorations!(ax), hidespines!(ax)
xlims!(-10, 7)
ylims!(-5, 7)

t = Node(1.0)
ecolorscaling = Node(0.5)
gpargs = gparguments(sol, t;
                     colortype=:abssteady,
                     offlinewidth=0.0,
                     ecolorscaling)
p = graphplot!(ax, network; gpargs...)
timelabel = @lift "time = " * repr(round($t, digits=2)) * " s"
fig[1,1] = Label(fig, timelabel, textsize=30,
                 halign=:left, valign=:top, tellwidth=false, tellheight=false)

largs = (;color=:black, linestyle=:dot, linewidth=3.0)

function recordsegment(fig, timeranges, i; line=false, fps=30)
    file = "rts_$i"
    tspan = timeranges[i]
    save(joinpath(DIR, "../videos", file*".png"), fig)
    record(fig, joinpath(DIR, "../videos", file*".mp4"), tspan; framerate=fps) do time
        t[] = time
        if line isa Int && time == timeranges[i][end]
            lines!(ax, getline(line); largs...)
        end
    end
end
l27 = lines!(ax, getline(27); largs...)

recordsegment(fig, timeranges, 1; line=29)

deleteat!(ax.scene.plots, findall(p->p isa Lines, ax.scene.plots))
division = lines!(ax, [(0,-3.0), (0.5,-5)]; largs...)

recordsegment(fig, timeranges, 2; line=83) # not importand
recordsegment(fig, timeranges, 3; line=37)
recordsegment(fig, timeranges, 4; fps=60)
